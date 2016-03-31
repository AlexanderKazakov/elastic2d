#include <lib/linal/special/Matrix22.hpp>

using namespace gcm;
using namespace gcm::linal;

real gcm::linal::
determinant(const Matrix22& m) {
	return determinant(m.xx, m.xy,
	                   m.yx, m.yy);
}


