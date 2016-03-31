#include <lib/linal/special/Matrix33.hpp>

using namespace gcm;
using namespace gcm::linal;

real gcm::linal::
determinant(const Matrix33& m) {
	return determinant(m.a11, m.a12, m.a13,
	                   m.a21, m.a22, m.a23,
	                   m.a31, m.a32, m.a33);
}


