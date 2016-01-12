#include "lib/util/areas/SphereArea.hpp"

using namespace gcm;


SphereArea::SphereArea(const real &radius, const linal::Vector3 &center) : radius(radius), center(center) {}

bool SphereArea::contains(const linal::Vector3 &coords) const {
	// TODO uncomment
	/*return (coords - center)*(coords - center) < radius*radius;*/
}
