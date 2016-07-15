#ifndef LIBGCM_UTILS_HPP
#define LIBGCM_UTILS_HPP

#include <cmath>
#include <algorithm>

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
//		assert_ne(t, 0); FIXME - return assert after skull calculation success
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
	
	
	/**
	 * Seed random generator to produce different values
	 */
	static void seedRand() {
		srand((unsigned int)time(0));
	}
	
	
	/**
	 * Produce pseudorandom uniformly distributed real number from min to max inclusive
	 * @note do not forget seedRand
	 */
	static real randomReal(const real min, const real max) {
		return ((max - min) * rand()) / RAND_MAX + min;
	}
	
	
	/**
	 * Check is the container has the value
	 */
	template<typename TContainer, typename TValue>
	static bool has(const TContainer& container, const TValue& value) {
		return std::find(container.begin(), container.end(), value) != container.end();
	}
	
	
	/**
	 * Return different number from {0, 1, 2}:
	 * 0,1 -> 2; 1,2 -> 0; 0,2 -> 1.
	 */
	static int other012(const int i, const int j) {
		assert_true(i != j && i >= 0 && i < 3 && j >= 0 && j < 3);
		for (int k = 0; k < 2; k++) {
			if (k != i && k != j) { return k; }
		}
		return 2;
	}
	
};


}

#endif /* LIBGCM_UTILS_HPP */
