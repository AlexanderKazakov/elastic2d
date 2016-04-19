#ifndef LIBGCM_UTILS_HPP
#define LIBGCM_UTILS_HPP

#include <cmath>

#include <lib/util/Types.hpp>
#include <lib/util/Assertion.hpp>

namespace gcm {
class Utils {
public:

	/**
	 * Signum function
	 * @throws Exception for zero argument
	 */
	template<typename T>
	static int sign(const T& t) {
		assert_ne(t, 0);
		return (t > 0) ? 1 : -1;
	}
	
	/**
	 * Approach to compare two real numbers with given tolerance.
	 * 
	 * For big numbers, it's same to 
	 * \f$   \abs{f1 - f2} < tolerance * \abs{f1 + f2} / 2   \f$.
	 * For small numbers, it's same to 
	 * \f$   \abs{f1 - f2} < tolerance^(3/2) / 2   \f$.
	 */
	static inline bool
	approximatelyEqual(const real f1, const real f2, 
	                   const real tolerance = EQUALITY_TOLERANCE) {
		
		real relativeError2 =  4 * (f1 - f2) * (f1 - f2) /
		                          ((f1 + f2) * (f1 + f2) + tolerance);
		
		return relativeError2 < tolerance * tolerance;
	}

};


}

#endif /* LIBGCM_UTILS_HPP */
