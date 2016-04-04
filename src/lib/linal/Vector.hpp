#ifndef LIBGCM_LINAL_VECTOR_HPP
#define LIBGCM_LINAL_VECTOR_HPP

#include <lib/linal/Matrix.hpp>

namespace gcm {
namespace linal {
/**
 * Generic vector - just a matrix with one column
 */
template<int TM, typename Container = DefaultMatrixContainer<TM, 1> >
using Vector = Matrix<TM, 1, Container>;

typedef Vector<1> Real1;
typedef Vector<2> Real2;
typedef Vector<3> Real3;


inline Real3 crossProduct(const Real3& v1, const Real3& v2) {
	return Real3({v1(1) * v2(2) - v1(2) * v2(1),
	              v1(2) * v2(0) - v1(0) * v2(2),
	              v1(0) * v2(1) - v1(1) * v2(0)});
}


inline Real2 perpendicularClockwise(const Real2& v) {
	return {v(1), -v(0)};
}


}
}

#endif // LIBGCM_LINAL_VECTOR_HPP
