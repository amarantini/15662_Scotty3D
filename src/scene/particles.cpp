
#include "particles.h"
#include<iostream>

bool Particles::Particle::update(const PT::Aggregate &scene, Vec3 const &gravity, const float radius, const float dt) {

	//A4T4: particle update

	// Compute the trajectory of this particle for the next dt seconds.

	float t=dt;
	do {
		// (1) Build a ray representing the particle's path as if it travelled at constant velocity.
		Ray ray(position, velocity);
		// (2) Intersect the ray with the scene and account for collisions. Be careful when placing
		// collision points using the particle radius. Move the particle to its next position.
		float ddt = 0;
		PT::Trace trace = scene.hit(ray);
		float speed = velocity.norm();
		float hit_time = trace.distance / speed;
		float rt = radius / speed; //time to travel the distance of radius
		if (trace.hit && !isnan(trace.distance) && trace.distance<=0){
			//collision at a negative time
			velocity = velocity - 2*dot(velocity,trace.normal)*trace.normal;
			ddt = 0;
		} else if(trace.hit && !isnan(hit_time) //&& hit_time<=t) {
				&& hit_time-rt>=0
				&& hit_time-rt<=t){
			//collision happens before end of timestep
			ddt = hit_time-rt;
			t -= ddt;
			position += velocity * ddt;
			velocity = velocity - 2*dot(velocity,trace.normal)*trace.normal;
		} else {
			// no collision OR
			// collision happens after end of timestep trace.distance>t
			position += velocity * t;
			ddt = t;
			t = 0;
		}
		// (3) Account for acceleration due to gravity after updating position.
		velocity += gravity * ddt;
		// (4) Repeat until the entire time step has been consumed.
	} while(t>0);
	
	// (5) Decrease the particle's age and return 'false' if it should be removed.
	age -= dt;
	return age > 0;
}

void Particles::advance(const PT::Aggregate& scene, const Mat4& to_world, float dt) {

	if(step_size < EPS_F) return;

	step_accum += dt;

	while(step_accum > step_size) {
		step(scene, to_world);
		step_accum -= step_size;
	}
}

void Particles::step(const PT::Aggregate& scene, const Mat4& to_world) {

	std::vector<Particle> next;
	next.reserve(particles.size());

	for(Particle& p : particles) {
		if(p.update(scene, gravity, radius, step_size)) {
			next.emplace_back(p);
		}
	}

	if(rate > 0.0f) {

		//helpful when emitting particles:
		float cos = std::cos(Radians(spread_angle) / 2.0f);

		//will emit particle i when i == time * rate
		//(i.e., will emit particle when time * rate hits an integer value.)
		//so need to figure out all integers in [current_step, current_step+1) * step_size * rate
		//compute the range:
		double begin_t = current_step * double(step_size) * double(rate);
		double end_t = (current_step + 1) * double(step_size) * double(rate);

		uint64_t begin_i = uint64_t(std::max(0.0, std::ceil(begin_t)));
		uint64_t end_i = uint64_t(std::max(0.0, std::ceil(end_t)));

		//iterate all integers in [begin, end):
		for (uint64_t i = begin_i; i < end_i; ++i) {
			//spawn particle 'i':

			float y = lerp(cos, 1.0f, rng.unit());
			float t = 2 * PI_F * rng.unit();
			float d = std::sqrt(1.0f - y * y);
			Vec3 dir = initial_velocity * Vec3(d * std::cos(t), y, d * std::sin(t));

			Particle p;
			p.position = to_world * Vec3(0.0f, 0.0f, 0.0f);
			p.velocity = to_world.rotate(dir);
			p.age = lifetime; //NOTE: could adjust lifetime based on index
			next.push_back(p);
		}
	}

	particles = std::move(next);
	current_step += 1;
}

void Particles::reset() {
	particles.clear();
	step_accum = 0.0f;
	current_step = 0;
	rng.seed(seed);
}

bool operator!=(const Particles& a, const Particles& b) {
	return a.gravity != b.gravity
	|| a.radius != b.radius
	|| a.initial_velocity != b.initial_velocity
	|| a.spread_angle != b.spread_angle
	|| a.lifetime != b.lifetime
	|| a.rate != b.rate
	|| a.step_size != b.step_size
	|| a.seed != b.seed;
}
