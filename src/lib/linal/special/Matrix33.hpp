#ifndef LIBGCM_LINAL_MATRIX33_HPP
#define LIBGCM_LINAL_MATRIX33_HPP

#include <lib/linal/Matrix.hpp>

namespace gcm {
namespace linal {
/**
 * Specialized value container for 3x3 matrix.
 */
class Matrix33Container {
public:
	static const int SIZE = 3 * 3;                 // size of storage in units of gcm::real
	union {
		real values[SIZE];
		struct {
			union {
				real xx;
				real a11;
			};
			union {
				real xy;
				real a12;
			};
			union {
				real xz;
				real a13;
			};
			union {
				real yx;
				real a21;
			};
			union {
				real yy;
				real a22;
			};
			union {
				real yz;
				real a23;
			};
			union {
				real zx;
				real a31;
			};
			union {
				real zy;
				real a32;
			};
			union {
				real zz;
				real a33;
			};
		};
	};
};


/**
 * Specialized 3x3 matrix implementation.
 */
typedef Matrix<3, 3, Matrix33Container> Matrix33;

/**
 * @return matrix determinant
 */
real determinant(const Matrix33& m);

/**
 * Arbitrary type 3x3 determinant
 */
template<
        typename T11, typename T12, typename T13,
        typename T21, typename T22, typename T23,
        typename T31, typename T32, typename T33
        >
auto determinant(const T11 &m11, const T12 &m12, const T13 &m13,
                 const T21 &m21, const T22 &m22, const T23 &m23,
                 const T31 &m31, const T32 &m32, const T33 &m33)
->decltype(m11 * m22 * m33) {
	return m11 * (m22 * m33 - m23 * m32) -
	       m12 * (m21 * m33 - m23 * m31) +
	       m13 * (m21 * m32 - m22 * m31);
}
}
}

#endif // LIBGCM_LINAL_MATRIX33_HPP
