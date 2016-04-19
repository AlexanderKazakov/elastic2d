#ifndef LIBLIBGCM_LINAL_ROTATION_MATRIX_HPP
#define LIBLIBGCM_LINAL_ROTATION_MATRIX_HPP

#include <lib/linal/Matrix.hpp>

namespace gcm {
namespace linal {
/**
 * Returns rotation matrix to perform rotation around X axis.
 *
 * @param angle Angle to rotate to.
 *
 * @return Rotation matrix.
 */
inline Matrix33 getXRotationMatrix(real angle) {
	return {1.0,  0.0,        0.0,
	        0.0,  cos(angle), sin(angle),
	        0.0, -sin(angle), cos(angle)};
}

/**
 * Returns rotation matrix to perform rotation around Y axis.
 *
 * @param angle Angle to rotate to.
 *
 * @return Rotation matrix.
 */
inline Matrix33 getYRotationMatrix(real angle) {
	return {cos(angle), 0.0, -sin(angle),
	        0.0,        1.0, 0.0,
	        sin(angle), 0.0, cos(angle)};
}

/**
 * Returns rotation matrix to perform rotation around Z axis.
 *
 * @param angle Angle to rotate to.
 *
 * @return Rotation matrix.
 */
inline Matrix33 getZRotationMatrix(real angle) {
	return {cos(angle),  sin(angle), 0.0,
	        -sin(angle), cos(angle), 0.0,
	        0.0,         0.0,        1.0};
}


}
}

#endif
