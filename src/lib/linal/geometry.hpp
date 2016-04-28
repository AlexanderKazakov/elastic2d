#ifndef LIBGCM_LINAL_GEOMETRY_HPP
#define LIBGCM_LINAL_GEOMETRY_HPP

#include <lib/linal/linearSystems.hpp>

namespace gcm {
namespace linal {

enum class POSITION {
	INSIDE,
	OUTSIDE,
	BORDER,
	FIRST_BORDER,
	SECOND_BORDER
};

/**
 * Computes barycentric coordinates of point q in triangle {a, b, c}
 * The order is {lambda_a, lambda_b, lambda_c}
 */
inline Real3 barycentricCoordinates(const Real2& a, const Real2& b, const Real2& c,
                                    const Real2& q) {
					
	/// @see https://en.wikipedia.org/wiki/Barycentric_coordinate_system 
	linal::Matrix22 T = {a(0) - c(0), b(0) - c(0),
						 a(1) - c(1), b(1) - c(1)};
	Real2 lambda = linal::solveLinearSystem(T, q - c);
	
	return {lambda(0), lambda(1), 1 - lambda(0) - lambda(1)};
}


/**
 * Find intersection of two lines (a1, a2) and (b1, b2).
 * Lines are infinite. Lines must not be collinear.
 */
inline Real2 linesIntersection(const Real2 a1, const Real2 a2,
                               const Real2 b1, const Real2 b2) {
	/// if use line representation \f$ \vec{r} = \vec{r_0} + t * \vec{tau} \f$,
	/// we can write SLE on parameters t
	Real2 tau1 = a1 - a2;
	Real2 tau2 = b1 - b2;
	
	Matrix22 A = {tau1(0), -tau2(0),
	              tau1(1), -tau2(1)};
	Real2 b = {(b1 - a1)(0),
	           (b1 - a1)(1)};
	Real2 t = solveLinearSystem(A, b);
	
	return a1 + t(0) * tau1;
}


/**
 * Is triangle on given points degenerate
 */
inline bool isDegenerate(const Real2 a, const Real2 b, const Real2 c) {
	Real2 l = a - b;
	Real2 m = c - b;
	return determinant(l(0), l(1), m(0), m(1)) == 0;
}


/** 
 * Is a orthogonal to b 
 */
template<int TM>
bool isPerpendicular(const Vector<TM>& a, const Vector<TM>& b) {
	return linal::dotProduct(a, b) == 0;
}


/**
 * Returns position of point q relative to the angle {a, b, c}.
 * b is the vertex of the triangle.
 * "inside" means that going from side b-a to side b-c counterclockwise,
 * we meet the point q.
 * "first border" means line a-b, "second border" means line c-b.
 * 
 *       c                          a          
 *        /                          /         
 *       /                          /          
 *      /                          /           
 *     /  inside                  /  outside    
 *    /___________               /___________ 
 *   b             a            b             c
 */
inline POSITION positionRelativeToAngle(
		const Real2& a, const Real2& b, const Real2& c, const Real2& q) {

	if (isDegenerate(a, b, c)) {
	
		if (isDegenerate(a, b, q)) {
			return (length(q - a) < length(q - c)) ? POSITION::FIRST_BORDER
			                                       : POSITION::SECOND_BORDER;
		} else {
			return (crossProduct(a - b, q - b) > 0) ? POSITION::INSIDE
			                                        : POSITION::OUTSIDE;
		}
		
	} else {
		Real3 lambda = barycentricCoordinates(a, b, c, q);
		
		if        (lambda(0) == 0 && lambda(2) > 0) {
			return POSITION::SECOND_BORDER;
			
		} else if (lambda(2) == 0 && lambda(0) > 0) {
			return POSITION::FIRST_BORDER;
			
		} else {
			bool isInsideTheLeastAngle = lambda(0) > 0 && lambda(2) > 0;			
			bool isAngleRight = crossProduct(a - b, c - b) > 0;
			
			if ( ( isAngleRight &&  isInsideTheLeastAngle) ||
			     (!isAngleRight && !isInsideTheLeastAngle) ) {
				return POSITION::INSIDE;
				
			} else {
				return POSITION::OUTSIDE;
			}
			
		}
	}
}


}
}

#endif // LIBGCM_LINAL_GEOMETRY_HPP
