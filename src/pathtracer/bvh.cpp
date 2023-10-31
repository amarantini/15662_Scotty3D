
#include "bvh.h"
#include "aggregate.h"
#include "instance.h"
#include "tri_mesh.h"

#include <stack>
#include <functional>
#include <iostream>

namespace PT {

struct BVHBuildData {
	BVHBuildData(size_t start, size_t range, size_t dst) : start(start), range(range), node(dst) {
	}
	size_t start; ///< start index into the primitive array
	size_t range; ///< range of index into the primitive array
	size_t node;  ///< address to update
};

struct SAHBucketData {
	BBox bb;          ///< bbox of all primitives
	size_t num_prims = 0; ///< number of primitives in the bucket
};

template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	//A3T3 - build a bvh

	// Keep these
    nodes.clear();
    primitives = std::move(prims);

    // Construct a BVH from the given vector of primitives and maximum leaf
    // size configuration.

	BBox bbox;
	for(size_t i=0; i<primitives.size(); i++){
		bbox.enclose(primitives[i].bbox());
	}

	// Declare lambda function first to enable recursion
	// std::function<void(int)> recursive_build;
	size_t idx_ = new_node(bbox,0,primitives.size(),0,0);

	auto recursive_build = [&](int idx, auto& recursive_build){
		
		size_t size = nodes[idx].size;
		size_t start = nodes[idx].start;
		// std::cout<<size<<","<<start<<"\n";
		float surface_area = nodes[idx].bbox.surface_area();
		if(size<=max_leaf_size){
			return;
		}

		int n_bucket = std::max(2, std::min(8,int(size/max_leaf_size)));
		SAHBucketData left_best, right_best;
		float best_cost = FLT_MAX;
		int best_axis = 0;
		for(int axis=0; axis<3; axis++){
			std::vector<SAHBucketData> bucket_list(n_bucket);
			sort(primitives.begin()+start, primitives.begin()+start+size, 
			[=](Primitive& a, Primitive& b){
				return a.bbox().center()[axis] < b.bbox().center()[axis];
			});

			float bucket_min = primitives[start].bbox().center()[axis];
			float bucket_max = primitives[start+size-1].bbox().center()[axis];
			float scale = float(n_bucket) / (bucket_max - bucket_min);

			for(size_t i=start; i<start+size; i++){
				BBox bbox = primitives[i].bbox();
				int idx = std::max(std::min(int((bbox.center()[axis]-bucket_min) * scale), n_bucket-1),0);
				bucket_list[idx].bb.enclose(bbox);
				bucket_list[idx].num_prims++;
			}
			
			SAHBucketData left;
			for(int j=1; j<n_bucket;j++){
				left.bb.enclose(bucket_list[j-1].bb);
				left.num_prims += bucket_list[j-1].num_prims;
				SAHBucketData right;
				for(int i=j; i<n_bucket; i++){
					right.bb.enclose(bucket_list[i].bb);
					right.num_prims += bucket_list[i].num_prims;
				}
				
				float cost = left.bb.surface_area() / surface_area * left.num_prims +
							right.bb.surface_area() / surface_area * right.num_prims;
				if(cost < best_cost){
					best_cost = cost;
					left_best = left;
					right_best = right;
					best_axis = axis;
				}
			}
		}
		

		sort(primitives.begin()+start, primitives.begin()+start+size, 
		[&](Primitive& a, Primitive& b){
			return a.bbox().center()[best_axis] < b.bbox().center()[best_axis];
		});
		nodes[idx].l = new_node(left_best.bb, start, left_best.num_prims,0,0);
		nodes[idx].r = new_node(right_best.bb, start+left_best.num_prims, right_best.num_prims,0,0);
		recursive_build(nodes[idx].l, recursive_build);
		recursive_build(nodes[idx].r, recursive_build);
	};
	
	recursive_build(idx_, recursive_build);

}

template<typename Primitive> Trace BVH<Primitive>::hit(const Ray& ray) const {
	//A3T3 - traverse your BVH

    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.

	//TODO: replace this code with a more efficient traversal:

	std::function<void(int,Trace&,Vec2)> find_closest_hit;
	find_closest_hit = [&](int idx, Trace& closest_hit, Vec2 times){

		if(nodes[idx].is_leaf()){
			for(size_t i=nodes[idx].start; i<nodes[idx].start+nodes[idx].size; i++) {
				Trace hit = primitives[i].hit(ray);
				closest_hit = Trace::min(closest_hit, hit);
			}
		} else {
			Vec2 times1 = times;
			Vec2 times2 = times;
			size_t l = nodes[idx].l;
			size_t r = nodes[idx].r;
			bool intersect1 = nodes[l].bbox.hit(ray, times1);
			bool intersect2 = nodes[r].bbox.hit(ray, times2);

			if(intersect1 && intersect2) {
				bool first = times1.x <= times2.x;
				if (first) {
					find_closest_hit(l, closest_hit, times1);
					if(!closest_hit.hit || times2.x <= closest_hit.distance) {
						find_closest_hit(r, closest_hit, times2);
					}
				} else {
					find_closest_hit(r, closest_hit, times2);
					if(!closest_hit.hit || times1.x <= closest_hit.distance) {
						find_closest_hit(l, closest_hit, times1);
					}
				}
			} else if (intersect1) {
				find_closest_hit(l, closest_hit, times1);
			} else if (intersect2) {
				find_closest_hit(r, closest_hit, times2);
			}
		}
	};

	Trace closest_hit;
	closest_hit.hit = false;
	closest_hit.distance = FLT_MAX;
	if(nodes.size()>0)
		find_closest_hit(root_idx, closest_hit, Vec2(0.0f,FLT_MAX));
    
    return closest_hit;
}

template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	build(std::move(prims), max_leaf_size);
}

template<typename Primitive> std::vector<Primitive> BVH<Primitive>::destructure() {
	nodes.clear();
	return std::move(primitives);
}

template<typename Primitive>
template<typename P>
typename std::enable_if<std::is_copy_assignable_v<P>, BVH<P>>::type BVH<Primitive>::copy() const {
	BVH<Primitive> ret;
	ret.nodes = nodes;
	ret.primitives = primitives;
	ret.root_idx = root_idx;
	return ret;
}

template<typename Primitive> Vec3 BVH<Primitive>::sample(RNG &rng, Vec3 from) const {
	if (primitives.empty()) return {};
	int32_t n = rng.integer(0, static_cast<int32_t>(primitives.size()));
	return primitives[n].sample(rng, from);
}

template<typename Primitive>
float BVH<Primitive>::pdf(Ray ray, const Mat4& T, const Mat4& iT) const {
	if (primitives.empty()) return 0.0f;
	float ret = 0.0f;
	for (auto& prim : primitives) ret += prim.pdf(ray, T, iT);
	return ret / primitives.size();
}

template<typename Primitive> void BVH<Primitive>::clear() {
	nodes.clear();
	primitives.clear();
}

template<typename Primitive> bool BVH<Primitive>::Node::is_leaf() const {
	// A node is a leaf if l == r, since all interior nodes must have distinct children
	return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
	Node n;
	n.bbox = box;
	n.start = start;
	n.size = size;
	n.l = l;
	n.r = r;
	nodes.push_back(n);
	return nodes.size() - 1;
}
 
template<typename Primitive> BBox BVH<Primitive>::bbox() const {
	if(nodes.empty()) return BBox{Vec3{0.0f}, Vec3{0.0f}};
	return nodes[root_idx].bbox;
}

template<typename Primitive> size_t BVH<Primitive>::n_primitives() const {
	return primitives.size();
}

template<typename Primitive>
uint32_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, uint32_t level,
                                   const Mat4& trans) const {

	std::stack<std::pair<size_t, uint32_t>> tstack;
	tstack.push({root_idx, 0u});
	uint32_t max_level = 0u;

	if (nodes.empty()) return max_level;

	while (!tstack.empty()) {

		auto [idx, lvl] = tstack.top();
		max_level = std::max(max_level, lvl);
		const Node& node = nodes[idx];
		tstack.pop();

		Spectrum color = lvl == level ? Spectrum(1.0f, 0.0f, 0.0f) : Spectrum(1.0f);
		GL::Lines& add = lvl == level ? active : lines;

		BBox box = node.bbox;
		box.transform(trans);
		Vec3 min = box.min, max = box.max;

		auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

		edge(min, Vec3{max.x, min.y, min.z});
		edge(min, Vec3{min.x, max.y, min.z});
		edge(min, Vec3{min.x, min.y, max.z});
		edge(max, Vec3{min.x, max.y, max.z});
		edge(max, Vec3{max.x, min.y, max.z});
		edge(max, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

		if (!node.is_leaf()) {
			tstack.push({node.l, lvl + 1});
			tstack.push({node.r, lvl + 1});
		} else {
			for (size_t i = node.start; i < node.start + node.size; i++) {
				uint32_t c = primitives[i].visualize(lines, active, level - lvl, trans);
				max_level = std::max(c + lvl, max_level);
			}
		}
	}
	return max_level;
}

template class BVH<Triangle>;
template class BVH<Instance>;
template class BVH<Aggregate>;
template BVH<Triangle> BVH<Triangle>::copy<Triangle>() const;

} // namespace PT
