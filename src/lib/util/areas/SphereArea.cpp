#include <lib/util/areas/SphereArea.hpp>

using namespace gcm;


SphereArea::SphereArea(const real& _radius, const linal::Vector3& _center) : radius(_radius), center(_center) {
	assert_gt(radius, 0.0);
}

bool SphereArea::contains(const linal::Vector3 &coords) const {
	return linal::length(coords - center) < radius;
}
