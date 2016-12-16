#ifndef LIBGCM_LINAL_GEOMETRY_HPP
#define LIBGCM_LINAL_GEOMETRY_HPP

#include <libgcm/linal/linearSystems.hpp>

namespace gcm {
namespace linal {


/**
 * Cross (aka vector) product
 */
inline Real3 crossProduct(const Real3& a, const Real3& b) {
	return Real3({a(1) * b(2) - a(2) * b(1),
	              a(2) * b(0) - a(0) * b(2),
	              a(0) * b(1) - a(1) * b(0)});
}


/**
 * Cross (aka vector) product for two vectors in XY-plane.
 * The result is equal to Z-component of cross product.
 * I.e, for right pair of vectors the result is positive,
 * for left pair - negative.
 */
inline real crossProduct(const Real2& a, const Real2& b) {
	return a(0) * b(1) - a(1) * b(0);
}


/** 
 * @return orthogonal to given vector of *the same length*
 * generated by clockwise rotation
 */
inline Real2 perpendicularClockwise(const Real2& v) {
	return {v(1), -v(0)};
}


/** 
 * @return orthogonal to given vector of *the same length*
 * clockwise perpendicular to the projection of given vector to XY-plane.
 * if given vector is orthogonal to XY-plane, i.e equal to (0, 0, v),
 * (v, 0, 0) is returned
 */
inline Real3 perpendicularClockwise(const Real3& v) {
	Real3 ans = {v(1), -v(0), 0};
	if (v(0) == 0 && v(1) == 0) {
		ans = {v(2), 0, 0};
	}
	return ans * length(v) / length(ans);
}


/**
 * Computes barycentric coordinates of point q in segment {a, b} in 1D.
 * The order is {lambda_a, lambda_b}
 */
inline Real2 barycentricCoordinates(
		const Real1& a, const Real1& b,
		const Real1& q) {
	real w = (q(0) - b(0)) / (a(0) - b(0));
	return {w, 1 - w};
}


/**
 * Computes barycentric coordinates of point q in segment {a, b} in 2D.
 * The order is {lambda_a, lambda_b}.
 * For correctness of the result point q must lie on the line {a, b}
 */
inline Real2 barycentricCoordinates(
		const Real2& a, const Real2& b,
		const Real2& q) {
	/// if q is on the line {a, b}, one of the lines of T is linearly dependent
	Matrix<3, 2> T = {   1,    1,
	                  a(0), b(0),
	                  a(1), b(1) };
	Real3 f = {1, q(0), q(1)};
	Real2 lambda = linearLeastSquares(T, f);
	return lambda / (lambda(0) + lambda(1));
}


/**
 * Computes barycentric coordinates of point q in segment {a, b} in 3D.
 * The order is {lambda_a, lambda_b}.
 * For correctness of the result point q must lie on the line {a, b}
 */
inline Real2 barycentricCoordinates(
		const Real3& a, const Real3& b,
		const Real3& q) {
	/// if q is on the line {a, b}, two of the lines of T is linearly dependent
	Matrix<4, 2> T = {   1,    1,
	                  a(0), b(0),
	                  a(1), b(1),
	                  a(2), b(2) };
	Real4 f = {1, q(0), q(1), q(2)};
	Real2 lambda = linearLeastSquares(T, f);
	return lambda / (lambda(0) + lambda(1));
}


/**
 * Computes barycentric coordinates of point q in triangle {a, b, c}
 * The order is {lambda_a, lambda_b, lambda_c}
 */
inline Real3 barycentricCoordinates(
		const Real2& a, const Real2& b, const Real2& c,
		const Real2& q) {
	/// @see https://en.wikipedia.org/wiki/Barycentric_coordinate_system 
	Matrix22 T = {a(0) - c(0), b(0) - c(0),
	              a(1) - c(1), b(1) - c(1)};
	Real2 lambda = solveLinearSystem(T, q - c);
	return {lambda(0), lambda(1), 1 - lambda(0) - lambda(1)};
}


/**
 * Computes barycentric coordinates of point q in triangle {a, b, c} in 3D
 * The order is {lambda_a, lambda_b, lambda_c}
 * For correctness of the result point q must lie in the flat {a, b, c}
 */
inline Real3 barycentricCoordinates(
		const Real3& a, const Real3& b, const Real3& c,
		const Real3& q) {
	/// if q is in the flat {a, b, c}, one of the lines of A is linearly dependent
	Matrix<4, 3> T = {   1,    1,    1,
	                  a(0), b(0), c(0), 
	                  a(1), b(1), c(1),
	                  a(2), b(2), c(2) };
	Real4 f = {1, q(0), q(1), q(2)};
	Real3 lambda = linearLeastSquares(T, f);
	return lambda / (lambda(0) + lambda(1) + lambda(2));
}


/**
 * Computes barycentric coordinates of point q in tetrahedron {a, b, c, d}
 * The order is {lambda_a, lambda_b, lambda_c, lambda_d}
 */
inline Real4 barycentricCoordinates(
		const Real3& a, const Real3& b, const Real3& c, const Real3& d,
		const Real3& q) {
	/// @see https://en.wikipedia.org/wiki/Barycentric_coordinate_system 
	Matrix33 T = {a(0) - d(0), b(0) - d(0), c(0) - d(0), 
	              a(1) - d(1), b(1) - d(1), c(1) - d(1),
	              a(2) - d(2), b(2) - d(2), c(2) - d(2),};
	Real3 lambda = solveLinearSystem(T, q - d);
	return {lambda(0), lambda(1), lambda(2), 1 - lambda(0) - lambda(1) - lambda(2)};
}


/**
 * Find intersection of two lines (a1, a2) and (b1, b2) in 2D.
 * Lines are infinite. Degenerate cases are not handled.
 */
inline Real2 linesIntersection(const Real2 a1, const Real2 a2,
                               const Real2 b1, const Real2 b2) {
	/// if use line representation \f$ \vec{r} = \vec{r_0} + t * \vec{tau} \f$,
	/// we can write SLE on parameters t
	Real2 tau1 = a1 - a2;
	Real2 tau2 = b1 - b2;
	
	Matrix22 A = {tau1(0), -tau2(0),
	              tau1(1), -tau2(1)};
	Real2 b = b1 - a1;
	Real2 t = solveLinearSystem(A, b);
	
	return a1 + t(0) * tau1;
}


/**
 * Find intersection of two lines (a1, a2) and (b1, b2) in 3D.
 * For correctness of the result, lines should intersect really (not checked)
 * Lines are infinite. Degenerate cases are not handled.
 */
inline Real3 linesIntersection(const Real3 a1, const Real3 a2,
                               const Real3 b1, const Real3 b2) {
	/// if use line representation \f$ \vec{r} = \vec{r_0} + t * \vec{tau} \f$,
	/// we can write SLE on parameters t
	Real3 tau1 = a1 - a2;
	Real3 tau2 = b1 - b2;
	
	Matrix<3, 2> A = {tau1(0), -tau2(0),
	                  tau1(1), -tau2(1),
	                  tau1(2), -tau2(2)};
	Real3 b = b1 - a1;
	Real2 t = linearLeastSquares(A, b);
	
	return a1 + t(0) * tau1;
}


/**
 * Find intersection of the flat specified by points f1, f2, f3
 * and the line specified by points l1, l2.
 * Line and flat are infinite. Degenerate cases are not handled.
 */
inline Real3 lineWithFlatIntersection(const Real3 f1, const Real3 f2, const Real3 f3,
		const Real3 l1, const Real3 l2) {
	/// if use line representation \f$ \vec{r} = \vec{l1} + k * \vec{tau} \f$,
	/// flat representation \f$ \vec{r} = \vec{f1} + m * \vec{p} + n * \vec{q} \f$,
	/// we can write SLE on parameters k, m, n
	Real3 tau = l2 - l1;
	Real3 p = f2 - f1;
	Real3 q = f3 - f1;
	
	Matrix33 A = {tau(0), -p(0), -q(0),
	              tau(1), -p(1), -q(1),
	              tau(2), -p(2), -q(2)};
	Real3 b = f1 - l1;
	Real3 params = solveLinearSystem(A, b);
	
	return l1 + params(0) * tau;
}


/**
 * Is triangle on given points degenerate with tolerance eps
 */
inline bool isDegenerate(
		const Real2 a, const Real2 b, const Real2 c, const real eps) {
	Real2 l = a - b;
	Real2 m = c - b;
	return std::fabs(determinant(l(0), l(1), m(0), m(1))) <= eps;
}


/**
 * Is tetrahedron on given points degenerate with tolerance eps
 */
inline bool isDegenerate(
		const Real3 a, const Real3 b, const Real3 c, const Real3 d, const real eps) {
	Real3 l = a - b;
	Real3 m = c - b;
	Real3 n = d - b;
	return std::fabs(determinant(
			l(0), l(1), l(2), m(0), m(1), m(2), n(0), n(1), n(2))) <= eps;
}


/**
 * Is triangle on given points degenerate with tolerance eps
 */
inline bool isDegenerate(
		const Real3 a, const Real3 b, const Real3 c, const real eps) {
	Real3 l = a - b;
	Real3 m = c - b;
	real lm = dotProduct(l, m);
	return std::fabs(lm * lm - dotProduct(l, l) * dotProduct(m, m)) <= eps;
}


/** 
 * @return oriented area of triangle {a, b, c};
 * Orientation is positive, if crossProduct(b - a, c - a) > 0.
 * In other words, if {a, b, c} is counterclockwise sequence 
 * looking from top of Z-axis.
 */
inline real orientedArea(const Real2 a, const Real2 b, const Real2 c) {
	return crossProduct(b - a, c - a) / 2;
}


/** @return area of triangle {a, b, c} */
inline real area(const Real2 a, const Real2 b, const Real2 c) {
	return fabs(orientedArea(a, b, c));
}


/** @return area of triangle {a, b, c} */
inline real area(const Real3 a, const Real3 b, const Real3 c) {
	return length(crossProduct(b - a, c - a)) / 2;
}


/** 
 * @return oriented volume of tetrahedron {a, b, c, d};
 * Orientation is positive, if dotProduct(d - a, crossProduct(b - a, c - a)) > 0.
 * In other words, if {a, b, c} is counterclockwise sequence looking from d.
 */
inline real orientedVolume(const Real3 a, const Real3 b, const Real3 c, const Real3 d) {
	Real3 ba = b - a;
	Real3 ca = c - a;
	Real3 da = d - a;
	return determinant( Matrix33({ba(0), ba(1), ba(2),
	                              ca(0), ca(1), ca(2),
	                              da(0), da(1), da(2)}) ) / 6;
}


/** @return volume of tetrahedron {a, b, c, d} */
inline real volume(const Real3 a, const Real3 b, const Real3 c, const Real3 d) {
	return fabs(orientedVolume(a, b, c, d));
}


/** 
 * Is a orthogonal to b 
 */
template<int TM>
bool isPerpendicular(const Vector<TM>& a, const Vector<TM>& b) {
	return linal::dotProduct(a, b) == 0;
}


/**
 * Does segment ab in 1D contain point q inside with tolerance eps
 */
inline bool segmentContains(Real1 a, Real1 b, Real1 q, real eps) {
	Real2 lambda = barycentricCoordinates(a, b, q);
	return lambda(0) >= -eps && lambda(1) >= -eps;
}


/**
 * Does segment ab in 2D contain point q inside with given tolerance.
 * Tolerance for degeneration test and for barycentric test can be different
 */
inline bool segmentContains(
		Real2 a, Real2 b, Real2 q, real eps, real degenerationEps) {
	if (!isDegenerate(a, b, q, degenerationEps)) { return false; }
	Real2 lambda = barycentricCoordinates(a, b, q);
	return lambda(0) >= -eps && lambda(1) >= -eps;
}


/**
 * Does segment ab in 3D contain point q inside with given tolerance.
 * Tolerance for degeneration test and for barycentric test can be different
 */
inline bool segmentContains(
		Real3 a, Real3 b, Real3 q, real eps, real degenerationEps) {
	if (!isDegenerate(a, b, q, degenerationEps)) { return false; }
	Real2 lambda = barycentricCoordinates(a, b, q);
	return lambda(0) >= -eps && lambda(1) >= -eps;
}


/**
 * Does triangle abc in 2D contain point q inside with tolerance eps
 */
inline bool triangleContains(Real2 a, Real2 b, Real2 c, Real2 q, real eps) {
	Real3 lambda = barycentricCoordinates(a, b, c, q);
	return lambda(0) >= -eps && lambda(1) >= -eps && lambda(2) >= -eps;
}


/**
 * Does triangle abc in 3D contain point q inside with given tolerance.
 * Tolerance for degeneration test and for barycentric test can be different
 */
inline bool triangleContains(
		Real3 a, Real3 b, Real3 c, Real3 q, real eps, real degenerationEps) {
	if (!isDegenerate(a, b, c, q, degenerationEps)) { return false; }
	Real3 lambda = barycentricCoordinates(a, b, c, q);
	return lambda(0) >= -eps && lambda(1) >= -eps && lambda(2) >= -eps;
}


/**
 * Does tetrahedron abcd contain point q inside with tolerance eps
 */
inline bool tetrahedronContains(
		const Real3& a, const Real3& b, const Real3& c, const Real3& d,
		const Real3& q, const real eps) {
	Real4 lambda = barycentricCoordinates(a, b, c, d, q);
	return lambda(0) >= -eps && lambda(1) >= -eps &&
	       lambda(2) >= -eps && lambda(3) >= -eps;
}


/**
 * Does angle b-a-c (a is in the middle) contain point q inside 
 * (i.e in its minimal sector) with tolerance eps.
 *           b/
 *  outside  /
 *          /
 *         /  inside
 *        /__________
 *       a           c
 */
inline bool angleContains(const Real2& a, const Real2& b, const Real2& c,
		const Real2& q, const real eps) {
	Real3 lambda = barycentricCoordinates(a, b, c, q);
	return lambda(0) <= 1 + eps && 
			lambda(1) >= -eps && lambda(2) >= -eps;
}


/**
 * Does solid angle a-{b-c-d} (a is in the center) contain point q inside
 * (i.e in its minimal sector) with tolerance eps.
 *           b/
 *  outside  /
 *          /
 *         /  inside
 *        /__________
 *       a\           c
 *         \
 *          \
 *          d
 */
inline bool solidAngleContains(
		const Real3& a, const Real3& b, const Real3& c, const Real3& d,
		const Real3& q, const real eps) {
	Real4 lambda = barycentricCoordinates(a, b, c, d, q);
	return lambda(0) <= 1 + eps && 
			lambda(1) >= -eps && lambda(2) >= -eps && lambda(3) >= -eps;
}


/**
 * Given with tetrahedron {opposite, a, b, c}, 
 * returns outer normal of tetrahedron's face {a, b, c}
 */
inline Real3 oppositeFaceNormal(const Real3& opposite, 
		const Real3& a, const Real3& b, const Real3& c) {

	Real3 ans = normalize(crossProduct(a - b, c - b));
	return (dotProduct(ans, a - opposite) > 0) ? ans : -ans;
}


/** Minimal height in triangle */
inline real minimalHeight(const Real2 a, const Real2 b, const Real2 c) {
	const real S = area(a, b, c);
	const real ab = length(a - b);
	const real ac = length(a - c);
	const real bc = length(b - c);
	return 2 * S / fmax(ab, fmax(ac, bc));
}


/** Minimal height in tetrahedron */
inline real minimalHeight(
		const Real3 a, const Real3 b, const Real3 c, const Real3 d) {
	const real V = volume(a, b, c, d);
	const real A = area(b, c, d);
	const real B = area(c, d, a);
	const real C = area(d, a, b);
	const real D = area(a, b, c);
	return 3 * V / fmax(A, fmax(B, fmax(C, D)));
}


/** 
 * Direction of the ray reflection from the surface
 * @param normal surface normal
 */
template<int TM>
Vector<TM> reflectionDirection(
		const Vector<TM> normal, const Vector<TM> initialDirection) {
	return normalize(initialDirection -
			2 * normal * dotProduct(initialDirection, normal));
}


}
}

#endif // LIBGCM_LINAL_GEOMETRY_HPP
