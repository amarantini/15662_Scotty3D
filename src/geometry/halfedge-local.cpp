
#include "halfedge.h"

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <iostream>

/******************************************************************
*********************** Local Operations **************************
******************************************************************/

/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it cannot perform an operation (i.e., because
    the resulting mesh does not have a valid representation).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementation, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/


/*
 * add_face: add a standalone face to the mesh
 *  sides: number of sides
 *  radius: distance from vertices to origin
 *
 * We provide this method as an example of how to make new halfedge mesh geometry.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::add_face(uint32_t sides, float radius) {
	//faces with fewer than three sides are invalid, so abort the operation:
	if (sides < 3) return std::nullopt;


	std::vector< VertexRef > face_vertices;
	//In order to make the first edge point in the +x direction, first vertex should
	// be at -90.0f - 0.5f * 360.0f / float(sides) degrees, so:
	float const start_angle = (-0.25f - 0.5f / float(sides)) * 2.0f * PI_F;
	for (uint32_t s = 0; s < sides; ++s) {
		float angle = float(s) / float(sides) * 2.0f * PI_F + start_angle;
		VertexRef v = emplace_vertex();
		v->position = radius * Vec3(std::cos(angle), std::sin(angle), 0.0f);
		face_vertices.emplace_back(v);
	}

	assert(face_vertices.size() == sides);

	//assemble the rest of the mesh parts:
	FaceRef face = emplace_face(false); //the face to return
	FaceRef boundary = emplace_face(true); //the boundary loop around the face

	std::vector< HalfedgeRef > face_halfedges; //will use later to set ->next pointers

	for (uint32_t s = 0; s < sides; ++s) {
		//will create elements for edge from a->b:
		VertexRef a = face_vertices[s];
		VertexRef b = face_vertices[(s+1)%sides];

		//h is the edge on face:
		HalfedgeRef h = emplace_halfedge();
		//t is the twin, lies on boundary:
		HalfedgeRef t = emplace_halfedge();
		//e is the edge corresponding to h,t:
		EdgeRef e = emplace_edge(false); //false: non-sharp

		//set element data to something reasonable:
		//(most ops will do this with interpolate_data(), but no data to interpolate here)
		h->corner_uv = a->position.xy() / (2.0f * radius) + 0.5f;
		h->corner_normal = Vec3(0.0f, 0.0f, 1.0f);
		t->corner_uv = b->position.xy() / (2.0f * radius) + 0.5f;
		t->corner_normal = Vec3(0.0f, 0.0f,-1.0f);

		//thing -> halfedge pointers:
		e->halfedge = h;
		a->halfedge = h;
		if (s == 0) face->halfedge = h;
		if (s + 1 == sides) boundary->halfedge = t;

		//halfedge -> thing pointers (except 'next' -- will set that later)
		h->twin = t;
		h->vertex = a;
		h->edge = e;
		h->face = face;

		t->twin = h;
		t->vertex = b;
		t->edge = e;
		t->face = boundary;

		face_halfedges.emplace_back(h);
	}

	assert(face_halfedges.size() == sides);

	for (uint32_t s = 0; s < sides; ++s) {
		face_halfedges[s]->next = face_halfedges[(s+1)%sides];
		face_halfedges[(s+1)%sides]->twin->next = face_halfedges[s]->twin;
	}

	return face;
}


/*
 * bisect_edge: split an edge without splitting the adjacent faces
 *  e: edge to split
 *
 * returns: added vertex
 *
 * We provide this as an example for how to implement local operations.
 * (and as a useful subroutine!)
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::bisect_edge(EdgeRef e) {
	// Phase 0: draw a picture
	//
	// before:
	//    ----h--->
	// v1 ----e--- v2
	//   <----t---
	//
	// after:
	//    --h->    --h2->
	// v1 --e-- vm --e2-- v2
	//    <-t2-    <--t--
	//

	// Phase 1: collect existing elements
	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;
	VertexRef v1 = h->vertex;
	VertexRef v2 = t->vertex;

	// Phase 2: Allocate new elements, set data
	VertexRef vm = emplace_vertex();
	vm->position = (v1->position + v2->position) / 2.0f;
	interpolate_data({v1, v2}, vm); //set bone_weights

	EdgeRef e2 = emplace_edge();
	e2->sharp = e->sharp; //copy sharpness flag

	HalfedgeRef h2 = emplace_halfedge();
	interpolate_data({h, h->next}, h2); //set corner_uv, corner_normal

	HalfedgeRef t2 = emplace_halfedge();
	interpolate_data({t, t->next}, t2); //set corner_uv, corner_normal

	// The following elements aren't necessary for the bisect_edge, but they are here to demonstrate phase 4
    FaceRef f_not_used = emplace_face();
    HalfedgeRef h_not_used = emplace_halfedge();

	// Phase 3: Reassign connectivity (careful about ordering so you don't overwrite values you may need later!)

	vm->halfedge = h2;

	e2->halfedge = h2;

	assert(e->halfedge == h); //unchanged

	//n.b. h remains on the same face so even if h->face->halfedge == h, no fixup needed (t, similarly)

	h2->twin = t;
	h2->next = h->next;
	h2->vertex = vm;
	h2->edge = e2;
	h2->face = h->face;

	t2->twin = h;
	t2->next = t->next;
	t2->vertex = vm;
	t2->edge = e;
	t2->face = t->face;
	
	h->twin = t2;
	h->next = h2;
	assert(h->vertex == v1); // unchanged
	assert(h->edge == e); // unchanged
	//h->face unchanged

	t->twin = h2;
	t->next = t2;
	assert(t->vertex == v2); // unchanged
	t->edge = e2;
	//t->face unchanged


	// Phase 4: Delete unused elements
    erase_face(f_not_used);
    erase_halfedge(h_not_used);

	// Phase 5: Return the correct iterator
	return vm;
}


/*
 * split_edge: split an edge and adjacent (non-boundary) faces
 *  e: edge to split
 *
 * returns: added vertex. vertex->halfedge should lie along e
 *
 * Note that when splitting the adjacent faces, the new edge
 * should connect to the vertex ccw from the ccw-most end of e
 * within the face.
 *
 * Do not split adjacent boundary faces.
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(EdgeRef e) {
	// A2L2 (REQUIRED): split_edge
	
	// collect
	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;

	FaceRef f1 = h->face;
	FaceRef f2 = t->face;

	VertexRef vm = bisect_edge(e).value();
	
	if(!f1->boundary){
		// create 
		FaceRef f_new = emplace_face(false);

		EdgeRef e_new = emplace_edge();
		HalfedgeRef h_new = emplace_halfedge();
		HalfedgeRef t_new = emplace_halfedge();
		h_new->twin = t_new;
		h_new->edge = e_new;

		t_new->twin = h_new;
		t_new->edge = e_new;

		e_new->halfedge = h_new;
		e_new->sharp = e->sharp;

		//collect
		HalfedgeRef hm1 = vm->halfedge;
		HalfedgeRef h_next_next = hm1->next->next;
		VertexRef v3 = h_next_next->vertex;

		// disconnect
		f1->halfedge = h_next_next;
		hm1->face = f_new;
		hm1->next->face = f_new;

		// conenct
		f_new->halfedge = hm1;

		h_new->vertex = v3;
		h_new->next = hm1;
		hm1->next->next = h_new;
		h_new->face = f_new;
		h_new->corner_normal = f_new->normal();

		t_new->vertex = vm;
		t_new->next = h_next_next;
		hm1->twin->next->twin->next = t_new;
		t_new->face = f1;
		t_new->corner_normal = f1->normal();
	}

	if(!f2->boundary) {
		// create 
		FaceRef f_new = emplace_face(false);
		EdgeRef e_new = emplace_edge();
		HalfedgeRef h_new = emplace_halfedge();
		HalfedgeRef t_new = emplace_halfedge();
		h_new->twin = t_new;
		h_new->edge = e_new;

		t_new->twin = h_new;
		t_new->edge = e_new;

		e_new->halfedge = h_new;
		e_new->sharp = e->sharp;

		//collect
		HalfedgeRef tm1 = vm->halfedge->twin->next;
		HalfedgeRef t_next_next = tm1->next->next;
		VertexRef v4 = t_next_next->vertex;

		// disconnect
		f2->halfedge = vm->halfedge->twin;
		tm1->face = f_new;
		tm1->next->face = f_new;

		// conenct
		f_new->halfedge = tm1;

		h_new->vertex = v4;
		h_new->next = tm1;
		tm1->next->next = h_new;
		h_new->corner_normal = f_new->normal();
		h_new->face = f_new;

		t_new->vertex = vm;
		t_new->next = t_next_next;
		vm->halfedge->twin->next = t_new;
		t_new->corner_normal = f2->normal();
		t_new->face = f2;
	}

    return vm;
}



/*
 * inset_vertex: divide a face into triangles by placing a vertex at f->center()
 *  f: the face to add the vertex to
 *
 * returns:
 *  std::nullopt if insetting a vertex would make mesh invalid
 *  the inset vertex otherwise
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::inset_vertex(FaceRef f) {
	// A2Lx4 (OPTIONAL): inset vertex
	if(f->boundary){
		return std::nullopt;
	}

	HalfedgeRef h_old = f->halfedge;
	std::vector<HalfedgeRef> old_hs;
	std::vector<HalfedgeCRef> old_hs_const;
	std::vector<VertexCRef> old_vs;
	std::vector<HalfedgeRef> new_hs;
	std::vector<FaceRef> new_fs = {f};
	int degree = f->degree();
	for(int i=0; i<degree; i++) {
		HalfedgeRef h = emplace_halfedge();
		HalfedgeRef t = emplace_halfedge();
		EdgeRef e = emplace_edge(false);
		if(i!=0) {
			FaceRef f = emplace_face(false);
			new_fs.push_back(f);
		}
		new_hs.push_back(h);

		h->twin = t;
		t->twin = h;
		h->edge = e;
		t->edge = e;
		e->halfedge = h;

		old_hs.push_back(h_old);
		old_hs_const.push_back(h_old);
		old_vs.push_back(h->vertex);
		h_old = h_old->next;
	}

	VertexRef v = emplace_vertex();
	v->halfedge = new_hs[0];

	for(int i=0; i<degree; i++) {
		new_fs[i]->halfedge = new_hs[i];

		new_hs[i]->vertex = v;
		new_hs[i]->next = old_hs[i];
		new_hs[i]->face = new_fs[i];

		int i_next = i==degree-1 ? 0: i+1;
		old_hs[i]->next = new_hs[i_next]->twin;
		old_hs[i]->face = new_fs[i];

		new_hs[i_next]->twin->next = new_hs[i];
		new_hs[i_next]->twin->face = new_fs[i];
		new_hs[i_next]->twin->vertex = old_hs[i]->twin->vertex;

		interpolate_data(old_hs_const, new_hs[i]);
		interpolate_data(old_hs_const, new_hs[i_next]->twin);
	}

	// interpolate_data(old_vs, v);
	v->position = v->neighborhood_center();
    return v;
}


/* [BEVEL NOTE] Note on the beveling process:

	Each of the bevel_vertex, bevel_edge, and extrude_face functions do not represent
	a full bevel/extrude operation. Instead, they should update the _connectivity_ of
	the mesh, _not_ the positions of newly created vertices. In fact, you should set
	the positions of new vertices to be exactly the same as wherever they "started from."

	When you click on a mesh element while in bevel mode, one of those three functions
	is called. But, because you may then adjust the distance/offset of the newly
	beveled face, we need another method of updating the positions of the new vertices.

	This is where bevel_positions and extrude_positions come in: these functions are
	called repeatedly as you move your mouse, the position of which determines the
	amount / shrink parameters. These functions are also passed an array of the original
	vertex positions, stored just after the bevel/extrude call, in order starting at
	face->halfedge->vertex, and the original element normal, computed just *before* the
	bevel/extrude call.

	Finally, note that the amount, extrude, and/or shrink parameters are not relative
	values -- you should compute a particular new position from them, not a delta to
	apply.
*/

/*
 * bevel_vertex: creates a face in place of a vertex
 *  v: the vertex to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(VertexRef v) {
	//A2Lx5 (OPTIONAL): Bevel Vertex
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in bevel_vertex_helper (A2Lx5h)

	(void)v;
    return std::nullopt;
}

/*
 * bevel_edge: creates a face in place of an edge
 *  e: the edge to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(EdgeRef e) {
	//A2Lx6 (OPTIONAL): Bevel Edge
	// Reminder: This function does not update the vertex positions.
	// remember to also fill in bevel_edge_helper (A2Lx6h)

	(void)e;
    return std::nullopt;
}

/*
 * extrude_face: creates a face inset into a face
 *  f: the face to inset
 *
 * returns: reference to the inner face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::extrude_face(FaceRef f) {
	//A2L4: Extrude Face
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in extrude_helper (A2L4h)

	// create
	uint32_t degree = f->degree();
	HalfedgeRef h = f->halfedge;

	std::vector<VertexRef> vs;
	std::vector<VertexRef> vs_outter;
	std::vector<FaceRef> fs;
	std::vector<HalfedgeRef> hs_inner;
	std::vector<HalfedgeRef> hs_connect;
	std::vector<HalfedgeRef> hs;

	HalfedgeRef htmp = h;
	for(uint32_t i=0; i<degree; i++) {
		VertexRef v = emplace_vertex();
		vs.push_back(v);

		FaceRef f_new = emplace_face(false);
		fs.push_back(f_new);

		// inner edge for f_new
		EdgeRef e1 = emplace_edge();
		HalfedgeRef h1 = emplace_halfedge();
		HalfedgeRef t1 = emplace_halfedge();
		e1->halfedge = h1;
		h1->twin = t1;
		t1->twin = h1;
		h1->edge = e1;
		t1->edge = e1;
		
		hs_inner.push_back(h1);

		// edge connecting v and old v
		EdgeRef e2 = emplace_edge();
		HalfedgeRef h2 = emplace_halfedge();
		HalfedgeRef t2 = emplace_halfedge();
		e2->halfedge = h2;
		h2->twin = t2;
		t2->twin = h2;
		h2->edge = e2;
		t2->edge = e2;
		
		hs_connect.push_back(h2);

		hs.push_back(htmp);
		vs_outter.push_back(htmp->vertex);
		htmp = htmp->next;
	}

	// connect
	
	f->halfedge = hs_inner[0];

	for(uint32_t i=0; i<degree; i++) {
		VertexRef v = vs[i];
		FaceRef f_new = fs[i];
		HalfedgeRef h1 = hs_inner[i];
		HalfedgeRef t1 = h1->twin;
		HalfedgeRef h2 = hs_connect[i+1==degree?0:i+1];
		HalfedgeRef t2 = hs_connect[i]->twin;
		HalfedgeRef htmp = hs[i];
		VertexRef v_o = vs_outter[i];
		v->position = v_o->position;

		v->halfedge = h1;
		
		h1->next = hs_inner[i+1==degree?0:i+1];
		h1->face = f;
		h1->corner_normal = htmp->corner_normal;
		h1->vertex = v;

		t1->face = f_new;
		t1->next = t2;
		t1->corner_normal = htmp->corner_normal;
		t1->vertex = vs[i+1==degree?0:i+1];

		t2->next = htmp;
		t2->face = f_new;
		t2->corner_normal = htmp->corner_normal;
		t2->vertex = v;
		
		htmp->next = h2;
		h2->next = t1;
		h2->face = f_new;
		h2->vertex = vs_outter[i+1==degree?0:i+1];
		h2->corner_normal = htmp->corner_normal;

		f_new->halfedge = htmp;
		htmp->face = f_new;
		v_o->halfedge = htmp;
	} 

	return f;
}

/*
 * flip_edge: rotate non-boundary edge ccw inside its containing faces
 *  e: edge to flip
 *
 * if e is a boundary edge, does nothing and returns std::nullopt
 * if flipping e would create an invalid mesh, does nothing and returns std::nullopt
 *
 * otherwise returns the edge, post-rotation
 *
 * does not create or destroy mesh elements.
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(EdgeRef e) {
	//A2L1: Flip Edge
	if(e->on_boundary())
    	return std::nullopt;

	// collect
	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;
	HalfedgeRef h_next = h->next;
	HalfedgeRef t_next = t->next;
	VertexRef v1 = h->next->vertex;
	VertexRef v2 = t->next->vertex;
	VertexRef v3 = h->next->next->vertex;
	VertexRef v4 = t->next->next->vertex;
	FaceRef f1 = h->face;
	FaceRef f2 = t->face;

	// edge case: a vertex only has 2 edges
	if(v1->degree()<=2 || v2->degree()<=2){
		return std::nullopt;
	}

	// disconnect
	v1->halfedge = h->next;
	v2->halfedge = t->next;
	f1->halfedge = h;
	f2->halfedge = t;

	// connect
	t->vertex = v3;
	h->vertex = v4;
	h->next = h_next->next;
	t->next = t_next->next;
	h_next->next = t;
	t_next->next = h;
	h_next->face = f2;
	t_next->face = f1;

	HalfedgeRef h_last = h;
	while(h_last->next != h){
		h_last = h_last->next;
	}
	h_last->next = t_next;

	HalfedgeRef t_last = t;
	while(t_last->next != t){
		t_last = t_last->next;
	}
	t_last->next = h_next;

	return e;
}


/*
 * make_boundary: add non-boundary face to boundary
 *  face: the face to make part of the boundary
 *
 * if face ends up adjacent to other boundary faces, merge them into face
 *
 * if resulting mesh would be invalid, does nothing and returns std::nullopt
 * otherwise returns face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::make_boundary(FaceRef face) {
	//A2Lx7: (OPTIONAL) make_boundary

	return std::nullopt; //TODO: actually write this code!
}

/*
 * dissolve_vertex: merge non-boundary faces adjacent to vertex, removing vertex
 *  v: vertex to merge around
 *
 * if merging would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_vertex(VertexRef v) {
	// A2Lx1 (OPTIONAL): Dissolve Vertex
	if(v->on_boundary() && v->degree()==2){
		return std::nullopt;
	}

	std::vector<HalfedgeRef> to_delete_h;
	std::vector<FaceRef> to_delete_f;
	std::vector<EdgeRef> to_delete_e;
	HalfedgeRef h = v->halfedge;
	HalfedgeRef h_b = halfedges.end(); //boundary halfedge
	HalfedgeRef h_b_twin_next = halfedges.end();
	FaceRef f = h->face;
	std::vector<HalfedgeRef> to_connect_h;
	do {
		if(h->face->boundary){
			h_b = h;
			h_b_twin_next = h->twin->next->twin->next->next;
		} else {
			to_delete_h.push_back(h);
			to_delete_h.push_back(h->twin);
			to_delete_e.push_back(h->edge);
			to_delete_f.push_back(h->face);
			to_connect_h.push_back(h->next);
			//disconnect
			h->next->face = f;
			h->twin->vertex->halfedge = h->next;
		}
	

		h = h->twin->next;
	} while (h != v->halfedge);

	//connect
	if(h_b != halfedges.end()){
		HalfedgeRef h_b_pprev = h_b;
		while(h_b_pprev->next->next!=h_b){
			h_b_pprev = h_b_pprev->next;
		}
		h_b->vertex = h_b_pprev->next->vertex;
		h_b_pprev->next = h_b;
		h_b->twin->face = f;
		h_b->twin->next = h_b_twin_next;
	}
		
	f->halfedge = to_connect_h.back();
	for(int i=to_connect_h.size()-1; i>=1; i--){
		if(h_b == halfedges.end() || to_connect_h[i]->next != h_b->twin)
			to_connect_h[i]->next = to_connect_h[i-1];
	}
	to_connect_h[0]->next = to_connect_h.back();

	Vec3 n = f->normal();
	for(int i=to_connect_h.size()-1; i>=1; i--){
		to_connect_h[i]->corner_normal = n;
	}

	// delete
	for(auto e: to_delete_e){
		erase_edge(e);
	}

	for(auto e: to_delete_h){
		erase_halfedge(e);
	}

	for(auto e: to_delete_f){
		if(e!=f)
			erase_face(e);
	}

	erase_vertex(v);

    return f;
}

/*
 * dissolve_edge: merge the two faces on either side of an edge
 *  e: the edge to dissolve
 *
 * merging a boundary and non-boundary face produces a boundary face.
 *
 * if the result of the merge would be an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_edge(EdgeRef e) {
	// A2Lx2 (OPTIONAL): dissolve_edge

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data

	// Edge case: dissolving an edge when there is only 1 face (+ 1 boundary face)
	if(e->on_boundary() && faces.size()==2)
    	return std::nullopt;

	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;
	FaceRef f1 = h->face;
	FaceRef f2 = t->face;

	// disconnect
	if(h->vertex->halfedge == h)
		h->vertex->halfedge = h->twin->next;

	if(t->vertex->halfedge == t)
		t->vertex->halfedge = t->twin->next;

	if(f2->boundary) {
		swap(h,t);
		swap(f1,f2);
	}

	// Merge f2 into f1
	do {
		t->face = f1;
		t = t->next;
	} while (t != h->twin);

	HalfedgeRef h_prev = h;
	do {
		h_prev = h_prev->next;
	} while (h_prev->next != h);

	HalfedgeRef t_prev = t;
	do {
		t_prev = t_prev->next;
	} while (t_prev->next != t);

	t_prev->next = h->next;
	h_prev->next = t->next;

	erase_edge(h->edge);
	erase_halfedge(h);
	erase_halfedge(t);
	erase_face(f2);

	return f1;
}

/* collapse_edge: collapse edge to a vertex at its middle
 *  e: the edge to collapse
 *
 * if collapsing the edge would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(EdgeRef e) {
	//A2L3: Collapse Edge

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)

	// collect
	HalfedgeRef h = e->halfedge;
    HalfedgeRef t = h->twin;

    VertexRef v1 = h->next->vertex;
    VertexRef v2 = t->next->vertex;
	FaceRef f1 = h->face;
    FaceRef f2 = t->face;
	uint32_t f1_degree = f1->degree();
	uint32_t f2_degree = f2->degree();


	auto joint_neighbor_edges_on_boundary = [=](VertexRef v1, VertexRef v2){
		// v1 and v2 has neighbor v3
		// if v1->v3 and v2->v3 on boundary, then return true
		HalfedgeRef h = v1->halfedge;
		do {
			if(h->next->next->vertex == v2) {
				return h->edge->on_boundary() && h->next->edge->on_boundary();
			}
			h = h->twin->next;
		} while (h != v1->halfedge);

		h = v2->halfedge;
		do {
			if(h->next->next->vertex == v1) {
				return h->edge->on_boundary() && h->next->edge->on_boundary();
			}
			h = h->twin->next;
		} while (h != v2->halfedge);
		return false; // no joint neighbor
	};

    if((!e->on_boundary() && 
		v1->on_boundary() && 
		v2->on_boundary() && 
		!joint_neighbor_edges_on_boundary(v1,v2))) {
        return std::nullopt;
	}
		

	auto merge_vertices = [&](){
		// move all edges of v2 to v1
		HalfedgeRef htmp = v2->halfedge;
		do {
			htmp->vertex = v1;
			htmp = htmp->twin->next;
		} while (htmp != v2->halfedge);
	};

	HalfedgeRef h_next = h->next;
	HalfedgeRef h_prev = h;
	while(h_prev->next != h){
		h_prev = h_prev->next;
	}
	HalfedgeRef t_next = t->next;
	HalfedgeRef t_prev = t;
	while(t_prev->next != t){
		t_prev = t_prev->next;
	}

	// disconnect

	v1->position = (v1->position + v2->position) / 2.0f;
	interpolate_data({v1, v2}, v1);
	std::set<HalfedgeRef> to_delete_he;
	to_delete_he.insert(h);
	to_delete_he.insert(t);

	merge_vertices();

	h_prev->next = h_next;
	t_prev->next = t_next;

	t_prev->face->halfedge = t_prev;

	auto merge_face = [=](FaceRef fa, std::set<HalfedgeRef>& to_delete_he){
		HalfedgeRef hm = fa->halfedge;
		HalfedgeRef hn = hm->next;
		FaceRef fb = hm->twin->face;

		// disconnect
		to_delete_he.insert(hm);
		to_delete_he.insert(hm->twin);
		while(to_delete_he.count(hm->twin->vertex->halfedge)){
			hm->twin->vertex->halfedge = hm->twin->vertex->halfedge->twin->next;
		}
		
		// if(hm->vertex->halfedge == hm)
		// 	hm->vertex->halfedge = hm->vertex->halfedge->twin->next;
		if(fb->halfedge == hm->twin)
			fb->halfedge = hm->twin->next;

		HalfedgeRef hm_prev = hm->next;
		if(hm->edge->on_boundary() && hm_prev->edge->on_boundary()){	
			HalfedgeRef tm = hm->twin;
			HalfedgeRef tm_prev = tm;
			do{
				tm_prev = tm_prev->next;
			} while(tm_prev->next != tm);

			to_delete_he.insert(hm_prev);
			to_delete_he.insert(hm_prev->twin);
			while(to_delete_he.count(hm_prev->vertex->halfedge)){
				hm_prev->vertex->halfedge = hm_prev->vertex->halfedge->twin->next;
			}

			tm_prev->next = tm->next->next;
			tm->face->halfedge = tm_prev;
			
			hm_prev->vertex->halfedge = hm_prev->twin->next;
			erase_edge(hm_prev->edge);
			erase_vertex(hm_prev->twin->vertex);
			erase_halfedge(hm_prev->twin);
			erase_halfedge(hm_prev);
			
		} else {
			// int count = 0;
			while(to_delete_he.count(hm->vertex->halfedge)){
				hm->vertex->halfedge = hm->vertex->halfedge->twin->next;
				// count++;
				// if(count>=20){
				// 	break;
				// }
			}

			// connect
			hn->next = hm->twin->next;
			hn->face = fb;
			HalfedgeRef hn_prev = fb->halfedge;
			do {
				hn_prev = hn_prev->next;
			} while(hn_prev->next != hm->twin);

			hn_prev->next = hn;
		}

		erase_face(fa);
		erase_edge(hm->edge);
		erase_halfedge(hm->twin);
		erase_halfedge(hm);
	};

	if(f1_degree>3 && f2_degree>3) {
		// no degenerating face
		while(to_delete_he.count(v1->halfedge))
			v1->halfedge = v1->halfedge->twin->next;

	} else if (f1_degree==3 && f2_degree==3){
		// remove both f1 and f2
		f1->halfedge = h_prev;
		f2->halfedge = t->next;
		
		merge_face(f1, to_delete_he);
		merge_face(f2, to_delete_he);
	} else if (f1_degree == 3) {
		// remove f1

		to_delete_he.insert(h_next);
		to_delete_he.insert(h_next->twin);
		to_delete_he.insert(h_prev);
		to_delete_he.insert(h_prev->twin);
		while(to_delete_he.count(v1->halfedge))
			v1->halfedge = v1->halfedge->twin->next;
		
		f1->halfedge = h_prev;
		merge_face(f1, to_delete_he);

	} else if (f2_degree == 3) {
		// remove f2

		to_delete_he.insert(h_next);
		to_delete_he.insert(h_next->twin);
		to_delete_he.insert(h_prev);
		to_delete_he.insert(h_prev->twin);
		while(to_delete_he.count(v1->halfedge))
			v1->halfedge = v1->halfedge->twin->next;
		
		f2->halfedge = t->next;
		merge_face(f2, to_delete_he);
	}

	erase_edge(e);
	erase_halfedge(t);
	erase_halfedge(h);
	erase_vertex(v2);

	return v1;
}

/*
 * collapse_face: collapse a face to a single vertex at its center
 *  f: the face to collapse
 *
 * if collapsing the face would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(FaceRef f) {
	//A2Lx3 (OPTIONAL): Collapse Face

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)
	if(f->boundary){
		return std::nullopt;
	}
	HalfedgeRef h = f->halfedge;

	// Edge case: if the face to collapse has >= 2 non-consequtive boundary edges, collapsing 
	// this face will result in hour glass shape
	auto has_non_consecutive_boundary_edges = [&](){
		int n = 0;
		do {
			if(h->edge->on_boundary())
				n++;
			do{
				h = h->next;
			} while(h->edge->on_boundary() && h!=f->halfedge);
		} while (h!=f->halfedge);
		return n>=2;
	};
	if(has_non_consecutive_boundary_edges()){
		return std::nullopt;
	}

	std::vector<HalfedgeRef> to_delete_h;
	std::vector<EdgeRef> to_delete_e;
	std::vector<VertexRef> to_delete_v;

	VertexRef v = h->vertex;
	if (v->halfedge->face == f || v->halfedge->twin->face==f) {
		v->halfedge = v->halfedge->twin->next;
	}

	do {
		to_delete_h.push_back(h);
		to_delete_h.push_back(h->twin);
		to_delete_e.push_back(h->edge);
		VertexRef vc = h->vertex;
		if(vc!=v){
			to_delete_v.push_back(h->vertex);
			HalfedgeRef hc = vc->halfedge;
			do {
				if(hc!=h){
					hc->vertex = v;
				}
				hc = hc->twin->next;
			} while (hc != vc->halfedge);
		}

		h = h->next;
	} while (h!=f->halfedge);

	do {
		HalfedgeRef t = h->twin;
		HalfedgeRef t_prev = t;
		interpolate_data({t_prev,t},t_prev);
		interpolate_data({t->next,t},t->next);
		do {
			t_prev = t_prev->next;
		} while (t_prev->next != t);
		t_prev->next = t->next;
		t->face->halfedge = t_prev;
		h = h->next;
	} while (h!=f->halfedge);
	
	// delete
	for(auto e: to_delete_e){
		erase_edge(e);
	}

	for(auto e: to_delete_h){
		erase_halfedge(e);
	}

	for(auto e: to_delete_v){
		if(e!=v)
			erase_vertex(e);
	}

	erase_face(f);

	// center v
	std::vector<VertexCRef> v_neighbors;
	h = v->halfedge;
	do {
		v_neighbors.push_back(h->twin->vertex);
		h = h->twin->next;
	} while(h != v->halfedge);

	v->position = v->neighborhood_center();
	interpolate_data(v_neighbors, v);

	return v;
}

/*
 * weld_edges: glue two boundary edges together to make one non-boundary edge
 *  e, e2: the edges to weld
 *
 * if welding the edges would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns e, updated to represent the newly-welded edge
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::weld_edges(EdgeRef e, EdgeRef e2) {
	//A2Lx8: Weld Edges

	//Reminder: use interpolate_data() to merge bone_weights data on vertices!

    return std::nullopt;
}



/*
 * bevel_positions: compute new positions for the vertices of a beveled vertex/edge
 *  face: the face that was created by the bevel operation
 *  start_positions: the starting positions of the vertices
 *     start_positions[i] is the starting position of face->halfedge(->next)^i
 *  direction: direction to bevel in (unit vector)
 *  distance: how far to bevel
 *
 * push each vertex from its starting position along its outgoing edge until it has
 *  moved distance `distance` in direction `direction`. If it runs out of edge to
 *  move along, you may choose to extrapolate, clamp the distance, or do something
 *  else reasonable.
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after bevel_vertex or bevel_edge.
 * (So you can assume the local topology is set up however your bevel_* functions do it.)
 *
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::bevel_positions(FaceRef face, std::vector<Vec3> const &start_positions, Vec3 direction, float distance) {
	//A2Lx5h / A2Lx6h (OPTIONAL): Bevel Positions Helper
	
	// The basic strategy here is to loop over the list of outgoing halfedges,
	// and use the preceding and next vertex position from the original mesh
	// (in the start_positions array) to compute an new vertex position.
	
}

/*
 * extrude_positions: compute new positions for the vertices of an extruded face
 *  face: the face that was created by the extrude operation
 *  move: how much to translate the face
 *  shrink: amount to linearly interpolate vertices in the face toward the face's centroid
 *    shrink of zero leaves the face where it is
 *    positive shrink makes the face smaller (at shrink of 1, face is a point)
 *    negative shrink makes the face larger
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after extrude_face.
 * (So you can assume the local topology is set up however your extrude_face function does it.)
 *
 * Using extrude face in the GUI will assume a shrink of 0 to only extrude the selected face
 * Using bevel face in the GUI will allow you to shrink and increase the size of the selected face
 * 
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::extrude_positions(FaceRef face, Vec3 move, float shrink) {
	//A2L4h: Extrude Positions Helper

	//General strategy:
	// use mesh navigation to get starting positions from the surrounding faces,
	// compute the centroid from these positions + use to shrink,
	// offset by move

	HalfedgeRef h = face->halfedge;

	Vec3 centroid = face->center();

	do {
		VertexRef v = h->vertex;
		if (shrink > 0.0f) {
			v->position = (v->position - centroid) * shrink + centroid;
		} else if (shrink < 0.0f) {
			v->position = (v->position - centroid) * 2.0f * (-1.0f) * shrink + centroid;
		}
		v->position += move;
		h = h->next;
	} while (h != face->halfedge);

	auto isNan = [=](Vec3 v){
		return v.x == NAN || v.y == NAN || v.z == NAN;
	};

	do {
		Vec3 n = face->normal();
		if(!isNan(n))
			h->corner_normal = n;
		HalfedgeRef t = h->twin;
		do {
			Vec3 n = face->normal();
			if(!isNan(n))
				t->corner_normal = n;
		} while(t != h->twin);
	} while(h != face->halfedge);
	
}

