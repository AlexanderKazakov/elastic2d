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
Matrix33 getXRotationMatrix(real angle);

/**
 * Returns rotation matrix to perform rotation around Y axis.
 *
 * @param angle Angle to rotate to.
 *
 * @return Rotation matrix.
 */
Matrix33 getYRotationMatrix(real angle);

/**
 * Returns rotation matrix to perform rotation around Z axis.
 *
 * @param angle Angle to rotate to.
 *
 * @return Rotation matrix.
 */
Matrix33 getZRotationMatrix(real angle);

}
}

#endif
