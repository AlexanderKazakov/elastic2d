#include <lib/util/areas/StraightBoundedCylinderArea.hpp>

using namespace gcm;


StraightBoundedCylinderArea::StraightBoundedCylinderArea(const real& _radius,
                                                         const Real3& _begin, const Real3& _end) :
		radius(_radius), begin(_begin), end(_end) {
	axis = linal::normalize(end - begin);
	assert_gt(radius, 0.0);
}

bool StraightBoundedCylinderArea::contains(const Real3 &coords) const {
	// check that the point is between caps
	if ( linal::dotProduct(coords - begin, axis) * linal::dotProduct(coords - end, axis) >= 0 ) return false;

	// check that the point is on the less than radius distance from the axis
	real projection = linal::dotProduct(coords - begin, axis);
	return linal::dotProduct(coords - begin, coords - begin) - projection*projection < radius*radius;
}

void StraightBoundedCylinderArea::move(const Real3& shift) {
	begin += shift; end += shift;
}