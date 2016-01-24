#include <lib/linal/special/Matrix22.hpp>

using namespace gcm;
using namespace gcm::linal;

real gcm::linal::determinant(const Matrix22 &m) {
	return m.xx * m.yy - m.xy * m.yx;
};

