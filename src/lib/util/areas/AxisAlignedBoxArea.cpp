#include <lib/util/areas/AxisAlignedBoxArea.hpp>

using namespace gcm;

AxisAlignedBoxArea::
AxisAlignedBoxArea(const Real3& _min, const Real3& _max) : min(_min), max(_max) {
	for (int i = 0; i < 3; i++) {
		assert_gt((max - min)(i), 0.0);
	}
}


bool AxisAlignedBoxArea::
contains(const Real3& coords) const {
	for (int i = 0; i < 3; i++) {
		if (coords(i) <= min(i) || coords(i) >= max(i)) {return false; }
	}
	return true;
}


void AxisAlignedBoxArea::
move(const Real3& shift) {
	min += shift; max += shift;
}


