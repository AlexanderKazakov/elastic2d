#include <lib/util/areas/StraightBoundedCylinderArea.hpp>

using namespace gcm;


StraightBoundedCylinderArea::StraightBoundedCylinderArea(const real &radius, const linal::Vector3 &begin, const linal::Vector3 &end) :
		radius(radius), begin(begin), end(end) {
	axis = linal::normalize(end - begin);
	assert_gt(radius, 0.0);
}

bool StraightBoundedCylinderArea::contains(const linal::Vector3 &coords) const {
	// check that the point is between caps
	if ( linal::dotProduct(coords - begin, axis) * linal::dotProduct(coords - end, axis) >= 0 ) return false;

	// check that the point is on the less than radius distance from the axis
	real projection = linal::dotProduct(coords - begin, axis);
	return linal::dotProduct(coords - begin, coords - begin) - projection*projection < radius*radius;
}
