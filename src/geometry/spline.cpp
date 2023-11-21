
#include "../geometry/spline.h"

template<typename T> T Spline<T>::at(float time) const {

	// A4T1b: Evaluate a Catumull-Rom spline

	// Given a time, find the nearest positions & tangent values
	// defined by the control point map.

	// Transform them for use with cubic_unit_spline

	// Be wary of edge cases! What if time is before the first knot,
	// before the second knot, etc...
	if(!any())
		return T();
	if(knots.size()==1){
		return knots.begin()->second;
	}
	if(time<=knots.begin()->first)
		return knots.begin()->second;
	if(time>=knots.rbegin()->first)
		return knots.rbegin()->second;

	auto k2_itr = knots.upper_bound(time); 
	auto k1_itr = std::prev(k2_itr);
	T k0,k3;
	T k1=k1_itr->second;
	T k2=k2_itr->second;
	float t0,t3;
	float t1 = k1_itr->first;
	float t2=k2_itr->first;
	if(k1_itr==knots.begin()){
		//don't have a knot "two to the left"
		k0 = k1 - (k2 - k1);
		t0 = t1 - (t2 - t1);
	} else {
		auto k0_itr = std::prev(k1_itr);
		k0 = k0_itr->second;
		t0 = k0_itr->first;
	}
	if(std::next(k2_itr)==knots.end()){
		//don't have a knot "two to the right"
		k3 = k2 + (k2 - k1);
		t3 = t2 + (t2 - t1);
	} else {
		auto k3_itr = std::next(k1_itr);
		k3 = k3_itr->second;
		t3 = k3_itr->first;
	}

	T m0 = (k2-k0) / (t2-t0);
	T m1 = (k3-k1) / (t3-t1);
	float t = (time - t1) / (t2 - t1);

	return cubic_unit_spline(t, k1, k2, m0, m1);
}

template<typename T>
T Spline<T>::cubic_unit_spline(float time, const T& position0, const T& position1,
                               const T& tangent0, const T& tangent1) {

	// A4T1a: Hermite Curve over the unit interval

	// Given time in [0,1] compute the cubic spline coefficients and use them to compute
	// the interpolated value at time 'time' based on the positions & tangents

	// Note that Spline is parameterized on type T, which allows us to create splines over
	// any type that supports the * and + operators.

	float t_2 = time * time;
	float t_3 = t_2 * time;

	float h00 = 2*t_3 - 3*t_2 + 1;
	float h10 = t_3 - 2*t_2 + time;
	float h01 = -2*t_3 + 3*t_2;
	float h11 = t_3 - t_2;


	return h00 * position0 + h10 * tangent0 + h01 * position1 + h11 * tangent1;
}

template class Spline<float>;
template class Spline<double>;
template class Spline<Vec4>;
template class Spline<Vec3>;
template class Spline<Vec2>;
template class Spline<Mat4>;
template class Spline<Spectrum>;
