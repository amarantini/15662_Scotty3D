// clang-format off
#include "pipeline.h"

#include <iostream>
#include <vector>

#include "../lib/log.h"
#include "../lib/mathlib.h"
#include "framebuffer.h"
#include "sample_pattern.h"
template<PrimitiveType primitive_type, class Program, uint32_t flags>
void Pipeline<primitive_type, Program, flags>::run(std::vector<Vertex> const& vertices,
                                                   typename Program::Parameters const& parameters,
                                                   Framebuffer* framebuffer_) {
	// Framebuffer must be non-null:
	assert(framebuffer_); 
	auto& framebuffer = *framebuffer_;

	// A1T7: sample loop
	// TODO: update this function to rasterize to *all* sample locations in the framebuffer.
	//  	 This will probably involve inserting a loop of the form:
	// 		 	std::vector< Vec3 > const &samples = framebuffer.sample_pattern.centers_and_weights;
	//      	for (uint32_t s = 0; s < samples.size(); ++s) { ... }
	//   	 around some subset of the code.
	// 		 You will also need to transform the input and output of the rasterize_* functions to
	// 	     account for the fact they deal with pixels centered at (0.5,0.5).
	std::vector< Vec3 > const &samples = framebuffer.sample_pattern.centers_and_weights;
	uint32_t sample_size = samples.size();
	for (uint32_t s = 0; s < sample_size; ++s) {

		std::vector<ShadedVertex> shaded_vertices;
		shaded_vertices.reserve(vertices.size());

	//--------------------------
	// shade vertices:
	for (auto const& v : vertices) {
		ShadedVertex sv;
		Program::shade_vertex(parameters, v.attributes, &sv.clip_position, &sv.attributes);
		shaded_vertices.emplace_back(sv);
	}

	//--------------------------
	// assemble + clip + homogeneous divide vertices:
	std::vector<ClippedVertex> clipped_vertices;

	// reserve some space to avoid reallocations later:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		// clipping lines can never produce more than one vertex per input vertex:
		clipped_vertices.reserve(shaded_vertices.size());
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		// clipping triangles can produce up to 8 vertices per input vertex:
		clipped_vertices.reserve(shaded_vertices.size() * 8);
	}
	// clang-format off

	//coefficients to map from clip coordinates to framebuffer (i.e., "viewport") coordinates:
	//x: [-1,1] -> [0,width]
	//y: [-1,1] -> [0,height]
	//z: [-1,1] -> [0,1] (OpenGL-style depth range)
	Vec3 const clip_to_fb_scale = Vec3{
		framebuffer.width / 2.0f,
		framebuffer.height / 2.0f,
		0.5f
	};
	Vec3 const clip_to_fb_offset = Vec3{
		0.5f * framebuffer.width,
		0.5f * framebuffer.height,
		0.5f
	};

	// helper used to put output of clipping functions into clipped_vertices:
	auto emit_vertex = [&](ShadedVertex const& sv) {
		ClippedVertex cv;
		float inv_w = 1.0f / sv.clip_position.w;
		cv.fb_position = clip_to_fb_scale * inv_w * sv.clip_position.xyz() + clip_to_fb_offset;
		cv.inv_w = inv_w;
		cv.attributes = sv.attributes;
		clipped_vertices.emplace_back(cv);
	};

	// actually do clipping:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		for (uint32_t i = 0; i + 1 < shaded_vertices.size(); i += 2) {
			clip_line(shaded_vertices[i], shaded_vertices[i + 1], emit_vertex);
		}
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		for (uint32_t i = 0; i + 2 < shaded_vertices.size(); i += 3) {
			clip_triangle(shaded_vertices[i], shaded_vertices[i + 1], shaded_vertices[i + 2], emit_vertex);
		}
	} else {
		static_assert(primitive_type == PrimitiveType::Lines, "Unsupported primitive type.");
	}

	//--------------------------
	// rasterize primitives:
	std::vector< Vec3 > const &samples = framebuffer.sample_pattern.centers_and_weights;
	for (uint32_t s = 0; s < samples.size(); ++s) {
	std::vector<Fragment> fragments;

	// helper used to put output of rasterization functions into fragments:
	auto emit_fragment = [&](Fragment const& f) { fragments.emplace_back(f); };

	// actually do rasterization:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		for (uint32_t i = 0; i + 1 < clipped_vertices.size(); i += 2) {
			auto cv1 = clipped_vertices[i];
			cv1.fb_position.x = cv1.fb_position.x - 0.5f + samples[s][0];
			cv1.fb_position.y = cv1.fb_position.y - 0.5f + samples[s][1];
			auto cv2 = clipped_vertices[i+1];
			cv2.fb_position.x = cv2.fb_position.x - 0.5f + samples[s][0];
			cv2.fb_position.y = cv2.fb_position.y - 0.5f + samples[s][1];
			rasterize_line(cv1, cv2, emit_fragment);
		}
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		for (uint32_t i = 0; i + 2 < clipped_vertices.size(); i += 3) {
			auto cv1 = clipped_vertices[i];
			cv1.fb_position.x = cv1.fb_position.x - 0.5f + samples[s][0];
			cv1.fb_position.y = cv1.fb_position.y - 0.5f + samples[s][1];
			auto cv2 = clipped_vertices[i+1];
			cv2.fb_position.x = cv2.fb_position.x - 0.5f + samples[s][0];
			cv2.fb_position.y = cv2.fb_position.y - 0.5f + samples[s][1];
			auto cv3 = clipped_vertices[i+2];
			cv3.fb_position.x = cv3.fb_position.x - 0.5f + samples[s][0];
			cv3.fb_position.y = cv3.fb_position.y - 0.5f + samples[s][1];
			rasterize_triangle(cv1, cv2, cv3, emit_fragment);
		}
	} else {
		static_assert(primitive_type == PrimitiveType::Lines, "Unsupported primitive type.");
	}

	//--------------------------
	// depth test + shade + blend fragments:
	uint32_t out_of_range = 0; // check if rasterization produced fragments outside framebuffer 
							   // (indicates something is wrong with clipping)
	for (auto const& f : fragments) {

		// fragment location (in pixels):
		int32_t x = (int32_t)std::floor(f.fb_position.x);
		int32_t y = (int32_t)std::floor(f.fb_position.y);

		// if clipping is working properly, this condition shouldn't be needed;
		// however, it prevents crashes while you are working on your clipping functions,
		// so we suggest leaving it in place:
		if (x < 0 || (uint32_t)x >= framebuffer.width || 
		    y < 0 || (uint32_t)y >= framebuffer.height) {
			++out_of_range;
			continue;
		}

		// local names that refer to destination sample in framebuffer:
		float& fb_depth = framebuffer.depth_at(x, y, s);
		Spectrum& fb_color = framebuffer.color_at(x, y, s);


		// depth test:
		if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Always) {
			// "Always" means the depth test always passes.
		} else if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Never) {
			// "Never" means the depth test never passes.
			continue; //discard this fragment
		} else if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Less) {
			// "Less" means the depth test passes when the new fragment has depth less than the stored depth.
			// A1T4: Depth_Less
			// TODO: implement depth test! We want to only emit fragments that have a depth less than the stored depth, hence "Depth_Less".
			if(fb_depth<=f.fb_position.z)
				continue;
		} else {
			static_assert((flags & PipelineMask_Depth) <= Pipeline_Depth_Always, "Unknown depth test flag.");
		}

		// if depth test passes, and depth writes aren't disabled, write depth to depth buffer:
		if constexpr (!(flags & Pipeline_DepthWriteDisableBit)) {
			fb_depth = f.fb_position.z;
		}

		// shade fragment:
		ShadedFragment sf;
		sf.fb_position = f.fb_position;
		Program::shade_fragment(parameters, f.attributes, f.derivatives, &sf.color, &sf.opacity);

		// write color to framebuffer if color writes aren't disabled:
		if constexpr (!(flags & Pipeline_ColorWriteDisableBit)) {
			// blend fragment:
			if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Replace) {
				fb_color = sf.color;
			} else if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Add) {
				// A1T4: Blend_Add
				// TODO: framebuffer color should have fragment color multiplied by fragment opacity added to it.
				fb_color += sf.color * sf.opacity; //<-- replace this line
			} else if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Over) {
				// A1T4: Blend_Over
				// TODO: set framebuffer color to the result of "over" blending (also called "alpha blending") the fragment color over the framebuffer color, using the fragment's opacity
				// 		 You may assume that the framebuffer color has its alpha premultiplied already, and you just want to compute the resulting composite color
				fb_color = fb_color * (1.0 - sf.opacity) + sf.color; //<-- replace this line
			} else {
				static_assert((flags & PipelineMask_Blend) <= Pipeline_Blend_Over, "Unknown blending flag.");
			}
		}
	}
	if (out_of_range > 0) {
		if constexpr (primitive_type == PrimitiveType::Lines) {
			warn("Produced %d fragments outside framebuffer; this indicates something is likely "
				"wrong with the clip_line function.",
				out_of_range);
		} else if constexpr (primitive_type == PrimitiveType::Triangles) {
			warn("Produced %d fragments outside framebuffer; this indicates something is likely "
				"wrong with the clip_triangle function.",
				out_of_range);
		}
	}
	}
	}
}

// -------------------------------------------------------------------------
// clipping functions

// helper to interpolate between vertices:
template<PrimitiveType p, class P, uint32_t F>
auto Pipeline<p, P, F>::lerp(ShadedVertex const& a, ShadedVertex const& b, float t) -> ShadedVertex {
	ShadedVertex ret;
	ret.clip_position = (b.clip_position - a.clip_position) * t + a.clip_position;
	for (uint32_t i = 0; i < ret.attributes.size(); ++i) {
		ret.attributes[i] = (b.attributes[i] - a.attributes[i]) * t + a.attributes[i];
	}
	return ret;
}

/*
 * clip_line - clip line to portion with -w <= x,y,z <= w, emit vertices of clipped line (if non-empty)
 *  	va, vb: endpoints of line
 *  	emit_vertex: call to produce truncated line
 *
 * If clipping shortens the line, attributes of the shortened line should respect the pipeline's interpolation mode.
 * 
 * If no portion of the line remains after clipping, emit_vertex will not be called.
 *
 * The clipped line should have the same direction as the full line.
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::clip_line(ShadedVertex const& va, ShadedVertex const& vb,
                                      std::function<void(ShadedVertex const&)> const& emit_vertex) {
	// Determine portion of line over which:
	// 		pt = (b-a) * t + a
	//  	-pt.w <= pt.x <= pt.w
	//  	-pt.w <= pt.y <= pt.w
	//  	-pt.w <= pt.z <= pt.w
	// ... as a range [min_t, max_t]:

	float min_t = 0.0f;
	float max_t = 1.0f;

	// want to set range of t for a bunch of equations like:
	//    a.x + t * ba.x <= a.w + t * ba.w
	// so here's a helper:
	auto clip_range = [&min_t, &max_t](float l, float dl, float r, float dr) {
		// restrict range such that:
		// l + t * dl <= r + t * dr
		// re-arranging:
		//  l - r <= t * (dr - dl)
		if (dr == dl) {
			// want: l - r <= 0
			if (l - r > 0.0f) {
				// works for none of range, so make range empty:
				min_t = 1.0f;
				max_t = 0.0f;
			}
		} else if (dr > dl) {
			// since dr - dl is positive:
			// want: (l - r) / (dr - dl) <= t
			min_t = std::max(min_t, (l - r) / (dr - dl));
		} else { // dr < dl
			// since dr - dl is negative:
			// want: (l - r) / (dr - dl) >= t
			max_t = std::min(max_t, (l - r) / (dr - dl));
		}
	};

	// local names for clip positions and their difference:
	Vec4 const& a = va.clip_position;
	Vec4 const& b = vb.clip_position;
	Vec4 const ba = b - a;

	// -a.w - t * ba.w <= a.x + t * ba.x <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.x, ba.x);
	clip_range(a.x, ba.x, a.w, ba.w);
	// -a.w - t * ba.w <= a.y + t * ba.y <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.y, ba.y);
	clip_range(a.y, ba.y, a.w, ba.w);
	// -a.w - t * ba.w <= a.z + t * ba.z <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.z, ba.z);
	clip_range(a.z, ba.z, a.w, ba.w);

	if (min_t < max_t) {
		if (min_t == 0.0f) {
			emit_vertex(va);
		} else {
			ShadedVertex out = lerp(va, vb, min_t);
			// don't interpolate attributes if in flat shading mode:
			if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
				out.attributes = va.attributes;
			}
			emit_vertex(out);
		}
		if (max_t == 1.0f) {
			emit_vertex(vb);
		} else {
			ShadedVertex out = lerp(va, vb, max_t);
			// don't interpolate attributes if in flat shading mode:
			if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
				out.attributes = va.attributes;
			}
			emit_vertex(out);
		}
	}
}

/*
 * clip_triangle - clip triangle to portion with -w <= x,y,z <= w, emit resulting shape as triangles (if non-empty)
 *  	va, vb, vc: vertices of triangle
 *  	emit_vertex: call to produce clipped triangles (three calls per triangle)
 *
 * If clipping truncates the triangle, attributes of the new vertices should respect the pipeline's interpolation mode.
 * 
 * If no portion of the triangle remains after clipping, emit_vertex will not be called.
 *
 * The clipped triangle(s) should have the same winding order as the full triangle.
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::clip_triangle(
	ShadedVertex const& va, ShadedVertex const& vb, ShadedVertex const& vc,
	std::function<void(ShadedVertex const&)> const& emit_vertex) {
	// A1EC: clip_triangle
	// TODO: correct code!
	emit_vertex(va);
	emit_vertex(vb);
	emit_vertex(vc);
}

// -------------------------------------------------------------------------
// rasterization functions

/*
 * rasterize_line:
 * calls emit_fragment( frag ) for every pixel "covered" by the line (va.fb_position.xy, vb.fb_position.xy).
 *
 *    a pixel (x,y) is "covered" by the line if it exits the inscribed diamond:
 * 
 *        (x+0.5,y+1)
 *        /        \
 *    (x,y+0.5)  (x+1,y+0.5)
 *        \        /
 *         (x+0.5,y)
 *
 *    to avoid ambiguity, we consider diamonds to contain their left and bottom points
 *    but not their top and right points. 
 * 
 * 	  since 45 degree lines breaks this rule, our rule in general is to rasterize the line as if its
 *    endpoints va and vb were at va + (e, e^2) and vb + (e, e^2) where no smaller nonzero e produces 
 *    a different rasterization result. 
 *    We will not explicitly check for 45 degree lines along the diamond edges (this will be extra credit),
 *    but you should be able to handle 45 degree lines in every other case (such as starting from pixel centers)
 *
 * for each such diamond, pass Fragment frag to emit_fragment, with:
 *  - frag.fb_position.xy set to the center (x+0.5,y+0.5)
 *  - frag.fb_position.z interpolated linearly between va.fb_position.z and vb.fb_position.z
 *  - frag.attributes set to va.attributes (line will only be used in Interp_Flat mode)
 *  - frag.derivatives set to all (0,0)
 *
 * when interpolating the depth (z) for the fragments, you may use any depth the line takes within the pixel
 * (i.e., you don't need to interpolate to, say, the closest point to the pixel center)
 *
 * If you wish to work in fixed point, check framebuffer.h for useful information about the framebuffer's dimensions.
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::rasterize_line(
	ClippedVertex const& va, ClippedVertex const& vb,
	std::function<void(Fragment const&)> const& emit_fragment) {
	if constexpr ((flags & PipelineMask_Interp) != Pipeline_Interp_Flat) {
		assert(0 && "rasterize_line should only be invoked in flat interpolation mode.");
	}
	// A1T2: rasterize_line
	float epsilon = 0.0001f;
	Vec2 a = va.fb_position.xy() + Vec2{epsilon,epsilon*epsilon};
	Vec2 b = vb.fb_position.xy() + Vec2{epsilon,epsilon*epsilon};
	float a_z = va.fb_position.z;
	float b_z = vb.fb_position.z;
	float d_x = b[0] - a[0];
	float d_y = b[1] - a[1];
	float d_z = vb.fb_position.z - va.fb_position.z;
	

	float x_a = va.fb_position.x;
	float y_a = va.fb_position.y;
	float x_b = vb.fb_position.x;
	float y_b = vb.fb_position.y;

	// Compute longer axis i and shorter axis j
	int i = 0, j = 1;
	if (std::abs(d_x) < std::abs(d_y)) {
		std::swap(i,j);
		std::swap(d_x,d_y);
		std::swap(x_a,y_a);
		std::swap(x_b,y_b);
	}

	float px_a = floor(x_a) + 0.5;
	float py_a = floor(y_a) + 0.5;
	float px_b = floor(x_b) + 0.5;
	float py_b = floor(y_b) + 0.5;

	// Starting coordinate has smaller value along longer axis
	if(a[i] > b[i]) {
		std::swap(a, b);
		d_x = -d_x;
		d_y = -d_y;
		d_z = -d_z;
		std::swap(a_z,b_z);
	}
	Vec2 ab = b-a;

	// Compute longer axis bound (exclude starting and ending point)
	float t1 = floor(a[i]);
	float t2 = floor(b[i]);

	// Shade start and end points
	/* - If the line starts inside a diamond and exits to continue on, 
	the starting pixel gets shaded. 
	- If the line starts at the outside of a diamond and continues on 
	without entering the said diamond, we do not shade the pixel. 
	- If the line ends inside a diamond, we do not shade the pixel. 
	- If the line ends outside a diamond it passed through, we shade the pixel.
	- A line to exit a diamond if 1) enters the diamond (either by starting in it 
	or passing through a point owned by the diamond) and 2) exits the diamond
	(by passing beyond any point on the diamond)*/

	auto pass_through = [&] (float px, float py) 
    {
		// line pass through diamond
        float w = (px - a[i])/d_x;
		float v = w * (d_y)+a[j];
		return std::abs(v-py)<0.5;
    };

	auto is_inside_diamond = [&] (float x, float y, float px, float py) 
    {
        return std::abs(y-py) + std::abs(x-px) < 0.5 ||
				(std::abs(px-x-0.5)<epsilon && std::abs(y-py)<epsilon)|| 
				(std::abs(py-y-0.5)<epsilon && std::abs(x-px)<epsilon);
    };
	
	bool is_a_inside_diamond = is_inside_diamond(x_a, y_a, px_a, py_a);
	bool is_b_inside_diamond = is_inside_diamond(x_b, y_b, px_b, py_b);
	
	bool shade_a = false, shade_b = false;
	if(std::abs(t1-t2)<epsilon){
		// a and b in the same pixel
		if(is_b_inside_diamond && is_a_inside_diamond)
			// a and b inside the same diamond
			return;
		else if (is_a_inside_diamond) {
			// a inside the diamond and b outside the diamond
			shade_a = true;
		} else if (!is_b_inside_diamond && !is_a_inside_diamond){
			// a outside the diamond and b outside the diamond
			shade_a = pass_through(px_a, py_a);
		}
	} else {
		if(is_b_inside_diamond)
			shade_a = is_a_inside_diamond;
		else {
			// a outside the diamond
			if(x_a <= px_a)
				shade_a = pass_through(px_a, py_a);
		}
		if(!is_b_inside_diamond) {
			// b outside the diamond
			if(x_b > px_b) // b pass through the diamond
				shade_b = pass_through(px_b, py_b);
		}
		
	}

	auto shade = [&](float px, float py){
		Vec2 pr = Vec2{px, py};
		float t = dot(pr - a,ab)/ab.norm();
		float z = t * (d_z)+a_z;
		Fragment frag;
		frag.fb_position = Vec3{px, py, z};;
		frag.attributes = va.attributes;
		frag.derivatives.fill(Vec2(0.0f, 0.0f));
		emit_fragment(frag);
	};

	if(shade_a) {
		if(i==1){
			std::swap(px_a,py_a);
		}
		shade(px_a,py_a);
	}

	if(shade_b) {
		if(i==1){
			std::swap(px_b,py_b);
		}
		shade(px_b,py_b);
	}
	
	
	// Ignore 2 endpoints
	t1 += 1.0;
	t2 -= 1.0;

	for(float u=t1; u<=t2; u+=1.0){
		float w = ((u+0.5) - a[i])/d_x;
		float v = w * (d_y)+a[j];
		float px = floor(u)+0.5, py = floor(v)+0.5;
		if(i==1){
			std::swap(px,py);
		}
		shade(px,py);
	}
	// The OpenGL specification section 3.5 may also come in handy.


}

/*
 * rasterize_triangle(a,b,c,emit) calls 'emit(frag)' at every location
 *  	(x+0.5,y+0.5) (where x,y are integers) covered by triangle (a,b,c).
 *
 * The emitted fragment should have:
 * - frag.fb_position.xy = (x+0.5, y+0.5)
 * - frag.fb_position.z = linearly interpolated fb_position.z from a,b,c (NOTE: does not depend on Interp mode!)
 * - frag.attributes = depends on Interp_* flag in flags:
 *   - if Interp_Flat: copy from va.attributes
 *   - if Interp_Smooth: interpolate as if (a,b,c) is a 2D triangle flat on the screen
 *   - if Interp_Correct: use perspective-correct interpolation
 * - frag.derivatives = derivatives w.r.t. fb_position.x and fb_position.y of the first frag.derivatives.size() attributes.
 *
 * Notes on derivatives:
 * 	The derivatives are partial derivatives w.r.t. screen locations. That is:
 *    derivatives[i].x = d/d(fb_position.x) attributes[i]
 *    derivatives[i].y = d/d(fb_position.y) attributes[i]
 *  You may compute these derivatives analytically or numerically.
 *
 *  See section 8.12.1 "Derivative Functions" of the GLSL 4.20 specification for some inspiration. (*HOWEVER*, the spec is solving a harder problem, and also nothing in the spec is binding on your implementation)
 *
 *  One approach is to rasterize blocks of four fragments and use forward and backward differences to compute derivatives.
 *  To assist you in this approach, keep in mind that the framebuffer size is *guaranteed* to be even. (see framebuffer.h)
 *
 * Notes on coverage:
 *  If two triangles are on opposite sides of the same edge, and a
 *  fragment center lies on that edge, rasterize_triangle should
 *  make sure that exactly one of the triangles emits that fragment.
 *  (Otherwise, speckles or cracks can appear in the final render.)
 * 
 *  For degenerate (co-linear) triangles, you may consider them to not be on any side of an edge.
 * 	Thus, even if two degnerate triangles share an edge that contains a fragment center, you don't need to emit it.
 *  You will not lose points for doing something reasonable when handling this case
 *
 *  This is pretty tricky to get exactly right!
 *
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::rasterize_triangle(
	ClippedVertex const& va, ClippedVertex const& vb, ClippedVertex const& vc,
	std::function<void(Fragment const&)> const& emit_fragment) {
	// NOTE: it is okay to restructure this function to allow these tasks to use the
	//  same code paths. Be aware, however, that all of them need to remain working!
	//  (e.g., if you break Flat while implementing Correct, you won't get points
	//   for Flat.)
	float epsilon = 0.00001f;
	Vec3 a = va.fb_position;
	Vec3 b = vb.fb_position;
	Vec3 c = vc.fb_position;

	if(std::abs(cross(b-a,c-a).norm())<epsilon) {
		return;
	}
	
	if(cross((c-a),(b-a)).z < 0){
		Vec3 tmp = c;
		c = b;
		b = tmp;
	}

	float left = floor(std::fmin(a[0],std::fmin(b[0],c[0])));
	float right = floor(std::fmax(a[0],std::fmax(b[0],c[0])));
	float bottom = floor(std::fmin(a[1],std::fmin(b[1],c[1])));
	float top = floor(std::fmax(a[1],std::fmax(b[1],c[1])));

	Vec3 top_edge, top_edge_v;
	Vec3 left_edge_1, left_edge_v_1;
	Vec3 left_edge_2, left_edge_v_2;
	bool has_top_edge = false;
	bool has_left_edge_1 = false;
	bool has_left_edge_2 = false;

	auto find_top_edge = [&]() {
		Vec3 ab = b-a;
		Vec3 bc = c-b;
		Vec3 ca = c-a;
		if(std::abs(ab[1])<epsilon && a[1]>c[1]){
			top_edge = ab;
			top_edge_v = a;
		} else if(std::abs(bc[1])<epsilon && b[1]>a[1]){
			top_edge = bc;
			top_edge_v = b;
		} else if(std::abs(ca[1])<epsilon && c[1]>b[1]){
			top_edge = ca;
			top_edge_v = c;
		} else {
			return;
		}
		has_top_edge = true;
	};

	auto find_left_edge = [&]() {
		Vec3 ab = b-a;
		Vec3 bc = c-b;
		Vec3 ca = a-c;
		if(ab[1]>0) {
			left_edge_1 = ab;
			left_edge_v_1 = a;
			has_left_edge_1 = true;
		}
		if(bc[1]>0) {
			if(!has_left_edge_1) {
				left_edge_1 = bc;
				left_edge_v_1 = b;
				has_left_edge_1 = true;
			} else {
				left_edge_2 = bc;
				left_edge_v_2 = b;
				has_left_edge_2 = true;
			}
			
		}
		if(ca[1]>0) {
			if(!has_left_edge_1) {
				left_edge_1 = ca;
				left_edge_v_1 = c;
				has_left_edge_1 = true;
			} else {
				left_edge_2 = ca;
				left_edge_v_2 = c;
				has_left_edge_2 = true;
			}
		}
		
	};

	find_top_edge();
	find_left_edge();

	auto on_edge = [=](Vec3 edge, Vec3 vert, Vec3 pt){
		return cross(pt-vert, edge)==Vec3();
	};	

	auto same_side = [=](Vec3& a, Vec3& b, Vec3& c, Vec3& pt){
		Vec3 ac_cross_ab = cross(c-a,b-a);
		Vec3 ac_cross_ap = cross(c-a,pt-a);
		return dot(ac_cross_ab,ac_cross_ap);
	};

	auto inside_triangle = [&](const float x, const float y, const float z)
	{   
		Vec3 pt = Vec3{x,y,z};
		float same_side_1 = same_side(a,b,c,pt);
		float same_side_2 = same_side(b,c,a,pt);
		float same_side_3 = same_side(c,a,b,pt);
		if(same_side_1>=0 && same_side_2>=0 && same_side_3>=0) {
			if (std::abs(same_side_1)<epsilon || std::abs(same_side_2)<epsilon || std::abs(same_side_3)<epsilon) {
				bool is_on_edge = (has_top_edge && on_edge(top_edge, top_edge_v, pt)) ||
					(has_left_edge_1 && on_edge(left_edge_1, left_edge_v_1, pt)) ||
					(has_left_edge_2 && on_edge(left_edge_2, left_edge_v_2, pt));
				return is_on_edge;
			}
			return true;
		}
		return false;
	};

	auto compute_barycentric = [=](float x, float y, Vec3 a, Vec3 b, Vec3 c)
	{
		float u = (x*(b[1] - c[1]) + (c[0] - b[0])*y + b[0]*c[1] - c[0]*b[1]) / (a[0]*(b[1] - c[1]) + (c[0] - b[0])*a[1] + b[0]*c[1] - c[0]*b[1]);
		float v = (x*(c[1] - a[1]) + (a[0] - c[0])*y + c[0]*a[1] - a[0]*c[1]) / (b[0]*(c[1] - a[1]) + (a[0] - c[0])*b[1] + c[0]*a[1] - a[0]*c[1]);
		float w = (x*(a[1] - b[1]) + (b[0] - a[0])*y + a[0]*b[1] - b[0]*a[1]) / (c[0]*(a[1] - b[1]) + (b[0] - a[0])*c[1] + a[0]*b[1] - b[0]*a[1]);
		return std::tuple{u,v,w};
	};

	auto interpolate = [=](float u, float v, float w, float x, float y, float z){
		return u*x + v*y + w*z;
	};

	auto pers_correct_interpolate = [=](float u,float v, float w, 
										float a_attr, float b_attr, float c_attr,
										float a_inv_w, float b_inv_w, float c_inv_w){
		float interpolated_attr_w_inv = interpolate(u,v,w,a_attr*a_inv_w,b_attr*b_inv_w,c_attr*c_inv_w);
		float interpolated_w_inv = interpolate(u,v,w,a_inv_w,b_inv_w,c_inv_w);
		return interpolated_attr_w_inv/interpolated_w_inv;
	};

	for(float i=left; i<=right; i++){
		for(float j=bottom; j<=top; j++) {
			float px = i+0.5;
			float py = j+0.5;
			auto[u, v, w] = compute_barycentric(px, py, va.fb_position, vb.fb_position, vc.fb_position);
			float pz = interpolate(u,v,w,va.fb_position.z,vb.fb_position.z,vc.fb_position.z);
			// float pz = pers_correct_interpolate(u,v,w,va.fb_position.z,vb.fb_position.z,vc.fb_position.z,va.inv_w,vb.inv_w,vc.inv_w);
			if(inside_triangle(px, py, pz)) {
				Fragment frag;
				frag.fb_position = Vec3{px, py, pz};
				
				if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
					// A1T3: flat triangles
					// TODO: rasterize triangle (see block comment above this function).

					// As a placeholder, here's code that draws some lines:
					//(remove this and replace it with a real solution)
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f)); 
				} else if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Smooth) {
					// A1T5: screen-space smooth triangles
					// TODO: rasterize triangle (see block comment above this function).

					// As a placeholder, here's code that calls the Flat interpolation version of the function:
					//(remove this and replace it with a real solution)
					//Forward difference - dx
					auto[u1, v1, w1] = compute_barycentric(px+1.0f, py, va.fb_position, vb.fb_position, vc.fb_position);
					//Forward difference - dy
					auto[u2, v2, w2] = compute_barycentric(px, py+1.0f, va.fb_position, vb.fb_position, vc.fb_position);
					for(long unsigned int i=0; i<va.attributes.size(); i++) {
						frag.attributes[i] = interpolate(u,v,w,va.attributes[i],vb.attributes[i],vc.attributes[i]);
					}
					for(long unsigned int i=0; i<2; i++) {
						frag.derivatives[i] = Vec2{
							interpolate(u1,v1,w1,va.attributes[i],vb.attributes[i],vc.attributes[i]) - frag.attributes[i],
							interpolate(u2,v2,w2,va.attributes[i],vb.attributes[i],vc.attributes[i]) - frag.attributes[i]
						};
					}
					
					
				} else if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Correct) {
					// A1T5: perspective correct triangles
					// TODO: rasterize triangle (block comment above this function).

					// As a placeholder, here's code that calls the Screen-space interpolation function:
					//(remove this and replace it with a real solution)
					//Forward difference - dx
					auto[u1, v1, w1] = compute_barycentric(px+1.0f, py, va.fb_position, vb.fb_position, vc.fb_position);
					//Forward difference - dy
					auto[u2, v2, w2] = compute_barycentric(px, py+1.0f, va.fb_position, vb.fb_position, vc.fb_position);
					for(long unsigned int i=0; i<va.attributes.size(); i++) {
						frag.attributes[i] = pers_correct_interpolate(u,v,w,va.attributes[i],vb.attributes[i],vc.attributes[i],va.inv_w,vb.inv_w,vc.inv_w);	
					}
					for(long unsigned int i=0; i<2; i++) {
						frag.derivatives[i] = Vec2{
							pers_correct_interpolate(u1,v1,w1,va.attributes[i],vb.attributes[i],vc.attributes[i],va.inv_w,vb.inv_w,vc.inv_w)
							- frag.attributes[i],
							pers_correct_interpolate(u2,v2,w2,va.attributes[i],vb.attributes[i],vc.attributes[i],va.inv_w,vb.inv_w,vc.inv_w)
							- frag.attributes[i]
						};
					}
					
				}
				emit_fragment(frag);
			}
		}
	}


	
}

//-------------------------------------------------------------------------
// compile instantiations for all programs and blending and testing types:

#include "programs.h"

template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Flat>;