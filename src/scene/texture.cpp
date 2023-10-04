
#include "texture.h"

#include <iostream>

namespace Textures {


Spectrum sample_nearest(HDR_Image const &image, Vec2 uv) {
	//clamp texture coordinates, convert to [0,w]x[0,h] pixel space:
	float x = image.w * std::clamp(uv.x, 0.0f, 1.0f);
	float y = image.h * std::clamp(uv.y, 0.0f, 1.0f);

	//the pixel with the nearest center is the pixel that contains (x,y):
	int32_t ix = int32_t(std::floor(x));
	int32_t iy = int32_t(std::floor(y));

	//texture coordinates of (1,1) map to (w,h), and need to be reduced:
	ix = std::min(ix, int32_t(image.w) - 1);
	iy = std::min(iy, int32_t(image.h) - 1);

	return image.at(ix, iy);
}

Spectrum sample_bilinear(HDR_Image const &image, Vec2 uv) {
	// A1T6: sample_bilinear
	//TODO: implement bilinear sampling strategy on texture 'image'
	float fx = int32_t(image.w) * std::clamp(uv.x, 0.0f, 1.0f);
	float fy = int32_t(image.h) * std::clamp(uv.y, 0.0f, 1.0f);
	int32_t w = int32_t(image.w);
	int32_t h = int32_t(image.h);
	
	int32_t x = std::max(std::min(static_cast<int32_t>(floor(fx-0.5f)),w-1),0);
	int32_t y = std::max(std::min(static_cast<int32_t>(std::floor(fy-0.5f)),h-1),0);
	// int32_t x = std::min(int32_t(std::floor(fx)),int32_t(image.w)-1);
	// int32_t y = std::min(int32_t(std::floor(fy)),int32_t(image.h)-1);
	int32_t x1,y1;

	// if(fx >= x+0.5f) {
	// 	x1 = std::min(x + 1, int32_t(image.w)-1);
	// } else {
		x1 = std::max(std::min(x + 1,w-1), 0);
	// }
	// if(fy>=y+0.5f){
	// 	y1 = std::min(y + 1, int32_t(image.h)-1);
	// } else {
		y1 = std::max(std::min(y + 1,h-1), 0);
	// }

	float dx = fx - 0.5f - float(x);
	float dy = fy - 0.5f - float(y);
	// float dx = std::abs(fx - std::floor(fx) - 0.5);
	// float dy = std::abs(fy - std::floor(fy) - 0.5);

	Spectrum t_xy = image.at(x,y);
	Spectrum t_x1y = image.at(x1,y);
	Spectrum t_xy1 = image.at(x,y1);
	Spectrum t_x1y1 = image.at(x1,y1);

	Spectrum tx = (1.0f - dx) * t_xy + dx * t_x1y;
	Spectrum ty = (1.0f - dx) * t_xy1 + dx * t_x1y1;

	return (1.0f - dy) * tx + dy * ty;
}


Spectrum sample_trilinear(HDR_Image const &base, std::vector< HDR_Image > const &levels, Vec2 uv, float lod) {
	// A1T6: sample_trilinear
	//TODO: implement trilinear sampling strategy on using mip-map 'levels'

	uint32_t d = std::min(std::max(0u,(uint32_t)floor(lod)), (uint32_t)(levels.size()));
	uint32_t d1 = std::min(std::max(0u,(uint32_t)floor(lod)+1), (uint32_t)(levels.size()));
	
	
	float dd = lod - floor(d);
	Spectrum td = sample_bilinear(d==0? base:levels[d-1], uv);
	Spectrum td1 = sample_bilinear(d1==0? base:levels[d1-1], uv);

	return (1.0f - dd) * td + dd * td1;
}

/*
 * generate_mipmap- generate mipmap levels from a base image.
 *  base: the base image
 *  levels: pointer to vector of levels to fill (must not be null)
 *
 * generates a stack of levels [1,n] of sizes w_i, h_i, where:
 *   w_i = max(1, floor(w_{i-1})/2)
 *   h_i = max(1, floor(h_{i-1})/2)
 *  with:
 *   w_0 = base.w
 *   h_0 = base.h
 *  and n is the smalles n such that w_n = h_n = 1
 *
 * each level should be calculated by downsampling a blurred version
 * of the previous level to remove high-frequency detail.
 *
 */
void generate_mipmap(HDR_Image const &base, std::vector< HDR_Image > *levels_) {
	assert(levels_);
	auto &levels = *levels_;


	{ // allocate sublevels sufficient to scale base image all the way to 1x1:
		int32_t num_levels = static_cast<int32_t>(std::log2(std::max(base.w, base.h)));
		assert(num_levels >= 0);

		levels.clear();
		levels.reserve(num_levels);

		uint32_t width = base.w;
		uint32_t height = base.h;
		for (int32_t i = 0; i < num_levels; ++i) {
			assert(!(width == 1 && height == 1)); //would have stopped before this if num_levels was computed correctly

			width = std::max(1u, width / 2u);
			height = std::max(1u, height / 2u);

			levels.emplace_back(width, height);
		}
		assert(width == 1 && height == 1);
		assert(levels.size() == uint32_t(num_levels));
	}

	//now fill in the levels using a helper:
	//downsample:
	// fill in dst to represent the low-frequency component of src
	auto downsample = [](HDR_Image const &src, HDR_Image &dst) {
		//dst is half the size of src in each dimension:
		assert(std::max(1u, src.w / 2u) == dst.w);
		assert(std::max(1u, src.h / 2u) == dst.h);

		// A1T6: generate
		//TODO: Write code to fill the levels of the mipmap hierarchy by downsampling
		for(int32_t i = 0; i < dst.w; i++) {
			for(int32_t j = 0; j < dst.h; j++) {
				Spectrum sum;
				for(int32_t x=2*i;x<2*(i+1) && x<src.w; x++){
					for(int32_t y=2*j;y<2*(j+1) && y<src.h; y++){
						sum += src.at(x,y);
					}
				}
				dst.at(i,j) = sum / 4.0f;
			}
		}
		//Be aware that the alignment of the samples in dst and src will be different depending on whether the image is even or odd.
	};

	std::cout << "Regenerating mipmap (" << levels.size() << " levels): [" << base.w << "x" << base.h << "]";
	std::cout.flush();
	for (uint32_t i = 0; i < levels.size(); ++i) {
		HDR_Image const &src = (i == 0 ? base : levels[i-1]);
		HDR_Image &dst = levels[i];
		std::cout << " -> [" << dst.w << "x" << dst.h << "]"; std::cout.flush();

		downsample(src, dst);
	}
	std::cout << std::endl;
	
}

Image::Image(Sampler sampler_, HDR_Image const &image_) {
	sampler = sampler_;
	image = image_.copy();
	update_mipmap();
}

Spectrum Image::evaluate(Vec2 uv, float lod) const {
	if (sampler == Sampler::nearest) {
		return sample_nearest(image, uv);
	} else if (sampler == Sampler::bilinear) {
		return sample_bilinear(image, uv);
	} else {
		return sample_trilinear(image, levels, uv, lod);
	}
}

void Image::update_mipmap() {
	if (sampler == Sampler::trilinear) {
		generate_mipmap(image, &levels);
	} else {
		levels.clear();
	}
}

GL::Tex2D Image::to_gl() const {
	return image.to_gl(1.0f);
}

void Image::make_valid() {
	update_mipmap();
}

Spectrum Constant::evaluate(Vec2 uv, float lod) const {
	return color * scale;
}

} // namespace Textures
bool operator!=(const Textures::Constant& a, const Textures::Constant& b) {
	return a.color != b.color || a.scale != b.scale;
}

bool operator!=(const Textures::Image& a, const Textures::Image& b) {
	return a.image != b.image;
}

bool operator!=(const Texture& a, const Texture& b) {
	if (a.texture.index() != b.texture.index()) return false;
	return std::visit(
		[&](const auto& data) { return data != std::get<std::decay_t<decltype(data)>>(b.texture); },
		a.texture);
}
