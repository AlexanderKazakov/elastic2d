#ifndef LIBGCM_LINAL_DETERMINANT_HPP
#define LIBGCM_LINAL_DETERMINANT_HPP

#include <lib/linal/Matrix.hpp>

namespace gcm {
namespace linal {

/**
 * Arbitrary type 2x2 determinant
 */
template<
        typename T11, typename T12,
        typename T21, typename T22
        >
inline auto determinant(const T11 &m11, const T12 &m12,
                        const T21 &m21, const T22 &m22)
->decltype(m11 * m22) {
	return m11 * m22 - m12 * m21;
}

inline real determinant(const Matrix22& m) {
	return determinant(m(0, 0), m(0, 1),
	                   m(1, 0), m(1, 1));
}


/**
 * Arbitrary type 3x3 determinant
 */
template<
        typename T11, typename T12, typename T13,
        typename T21, typename T22, typename T23,
        typename T31, typename T32, typename T33
        >
inline auto determinant(const T11 &m11, const T12 &m12, const T13 &m13,
                        const T21 &m21, const T22 &m22, const T23 &m23,
                        const T31 &m31, const T32 &m32, const T33 &m33)
->decltype(m11 * m22 * m33) {
	return m11 * (m22 * m33 - m23 * m32) -
	       m12 * (m21 * m33 - m23 * m31) +
	       m13 * (m21 * m32 - m22 * m31);
}

inline real determinant(const Matrix33& m) {
	return determinant(m(0, 0), m(0, 1), m(0, 2),
	                   m(1, 0), m(1, 1), m(1, 2),
	                   m(2, 0), m(2, 1), m(2, 2));
}


}
}


#endif // LIBGCM_LINAL_DETERMINANT_HPP
