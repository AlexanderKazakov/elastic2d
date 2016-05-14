#ifndef LIBGCM_LINAL_BASIS_HPP
#define LIBGCM_LINAL_BASIS_HPP

#include <lib/linal/functions.hpp>
#include <lib/linal/geometry.hpp>

namespace gcm {
namespace linal {

/** 
 * @name Local basis creation
 * Create local orthogonal basis for given vector n of unit length.
 * @warning vector n MUST BE of unit length
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
 * @return orthogonal transfer matrix from created basis to global {X, Y, Z} basis,
 * i.e matrix with vectors of created basis in columns.
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
 ///@}


}
}

#endif // LIBGCM_LINAL_BASIS_HPP
