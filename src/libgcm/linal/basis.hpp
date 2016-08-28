#ifndef LIBGCM_LINAL_BASIS_HPP
#define LIBGCM_LINAL_BASIS_HPP

#include <libgcm/linal/functions.hpp>
#include <libgcm/linal/geometry.hpp>
#include <libgcm/linal/special/RotationMatrix.hpp>

namespace gcm {
namespace linal {

/** 
 * @name Local basis creation
 * Create local orthogonal basis for given vector n of unit length.
 * @note vector n MUST BE of unit length
 * 
 * In created basis, given vector n is always at last position
 * and vector at the first position tau_1 is always in XY-plane.
 * {tau_1, tau_2, n} is right-hand orthogonal triple of vectors.
 * 
 * In 3D, it's impossible to go around the whole sphere with local basis,
 * continuos at every point. In current realization, the local basis is 
 * continuos while going around the sphere slice by XY-plane, but it has 
 * discontinuities while going through the points with normal(0, 0, n).
 * 
 * 
 *          Z |      tau_2
 *            |    / 
 *  n         |   / 
 *    \__     |  /
 *       \__  | /
 *          \ |/            Y
 *            .-\-----------
 *           / \/
 *          /   \
 *         /     \
 *        /       \ 
 *     X /          tau_1
 * 
 * 
 * Functions createLocalBasis return orthogonal transfer matrix 
 * from created basis to global {X, Y, Z} basis,
 * i.e matrix with vectors of created basis in columns.
 * 
 * Functions createLocalBasisTranspose return orthogonal transfer matrix 
 * from global {X, Y, Z} basis to created basis,
 * i.e matrix with vectors of created basis in strings.
 */
 ///@{
inline Matrix11 createLocalBasis(const Real1& n) {
	return Matrix11({n(0)});
}

inline Matrix22 createLocalBasis(const Real2& n) {
	const Real2 tau = perpendicularClockwise(n);
	return Matrix22({tau(0), n(0),
	                 tau(1), n(1)});
}

inline Matrix33 createLocalBasis(const Real3& n) {
	const Real3 tau_1 = normalize(perpendicularClockwise(n));
	const Real3 tau_2 = crossProduct(n, tau_1);
	return Matrix33({tau_1(0), tau_2(0), n(0),
	                 tau_1(1), tau_2(1), n(1),
	                 tau_1(2), tau_2(2), n(2)});
}


inline Matrix11 createLocalBasisTranspose(const Real1& n) {
	return Matrix11({n(0)});
}

inline Matrix22 createLocalBasisTranspose(const Real2& n) {
	const Real2 tau = perpendicularClockwise(n);
	return Matrix22({tau(0), tau(1),
	                 n(0),   n(1)});
}

inline Matrix33 createLocalBasisTranspose(const Real3& n) {
	const Real3 tau_1 = normalize(perpendicularClockwise(n));
	const Real3 tau_2 = crossProduct(n, tau_1);
	return Matrix33({tau_1(0), tau_1(1), tau_1(2),
	                 tau_2(0), tau_2(1), tau_2(2),
	                     n(0),     n(1),     n(2)});
}
 ///@}

/// @name Random basis creation
 /// @{
inline Matrix11 randomBasis(const Matrix11&) {
	return Matrix11({Utils::randomReal(-1, 1) > 0 ? 1.0 : -1.0});
}

inline Matrix22 randomBasis(const Matrix22&) {
	real phi = Utils::randomReal(-M_PI, M_PI);
	return createLocalBasis(normalize(Real2({cos(phi), sin(phi)})));
}

inline Matrix33 randomBasis(const Matrix33&) {
	real phi = Utils::randomReal(-M_PI, M_PI);
	real teta = Utils::randomReal(-M_PI, M_PI);
	real khi = Utils::randomReal(-M_PI, M_PI);
	
	Matrix33 ans = getZRotationMatrix(khi) * 
			getYRotationMatrix(teta) * getXRotationMatrix(phi);
	
	for (int i = 0; i < 3; i++) {
		ans.setColumn(i, normalize(ans.getColumn(i)));
	}
	return ans;
}
 ///@}


}
}

#endif // LIBGCM_LINAL_BASIS_HPP
