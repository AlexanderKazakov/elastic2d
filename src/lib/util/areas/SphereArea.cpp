#include <lib/util/areas/SphereArea.hpp>

using namespace gcm;


SphereArea::SphereArea(const real& _radius, const Real3& _center) : radius(_radius), center(_center) {
	assert_gt(radius, 0.0);
}

bool SphereArea::contains(const Real3 &coords) const {
	return linal::length(coords - center) < radius;
}

void SphereArea::move(const Real3& shift) {
	center += shift;
}