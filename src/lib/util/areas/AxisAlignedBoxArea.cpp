#include "lib/util/areas/AxisAlignedBoxArea.hpp"

using namespace gcm;

AxisAlignedBoxArea::AxisAlignedBoxArea(const linal::Vector3 &min, const linal::Vector3 &max) : min(min), max(max) {}

bool AxisAlignedBoxArea::contains(const linal::Vector3 &coords) const {
	for (int i = 0; i < 3; i++) {
		if (coords(i) <= min(i) || coords(i) >= max(i)) return false;
	}
	return true;
}
