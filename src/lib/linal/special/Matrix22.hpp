#ifndef LIBGCM_LINAL_MATRIX22_HPP
#define LIBGCM_LINAL_MATRIX22_HPP

#include <lib/linal/Matrix.hpp>

namespace gcm {
namespace linal {
/**
 * Specialized value container for 2x2 matrix.
 */
struct Matrix22Container {
	static const int SIZE = 2 * 2; ///< size of storage in units of gcm::real
	union {
		real values[SIZE];
		struct {
			real xx;
			real xy;
			real yx;
			real yy;
		};
	};
};


/**
 * Specialized 2x2 matrix implementation.
 */
typedef Matrix<2, 2, Matrix22Container> Matrix22;

/**
 * @return matrix determinant
 */
real determinant(const Matrix22& m);

/**
 * Arbitrary type determinant
 */
template<
        typename T11, typename T12,
        typename T21, typename T22
        >
auto determinant(const T11 &m11, const T12 &m12,
                 const T21 &m21, const T22 &m22)
->decltype(m11 * m22) {
	return m11 * m22 - m12 * m21;
};
};
};

#endif // LIBGCM_LINAL_MATRIX22_HPP
