#include "lib/util/areas/StraightBoundedCylinderArea.hpp"

using namespace gcm;


StraightBoundedCylinderArea::StraightBoundedCylinderArea(const real &radius, const linal::Vector3 &begin, const linal::Vector3 &end) :
		radius(radius), begin(begin), end(end) {
	axis = linal::normalize(end - begin);
}

bool StraightBoundedCylinderArea::contains(const linal::Vector3 &coords) const {
	// TODO uncomment
/*	// check that the point is between caps
	if ( ((coords - begin) * axis) * ((coords - end) * axis) > 0 ) return false;

	// check that the point is on the less than radius distance from the axis
	real projection = (coords - begin) * axis;
	return (coords - begin)*(coords - begin) - projection*projection < radius*radius;*/
}
