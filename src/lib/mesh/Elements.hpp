#ifndef LIBGCM_ELEMENTS_HPP
#define LIBGCM_ELEMENTS_HPP

#include <lib/linal/linal.hpp>
#include <set>

namespace gcm {
namespace elements {
/**
 * Convex geometric primitives by sets of their vertices
 * @tparam Point representation of point
 * @tparam TN number of points
 */
template<typename Point, int TN>
struct Element {

	static const int N = TN;
	Point p[N];
	bool inner = false;
	
	bool operator==(const Element& other) const {
		if (this->inner != other.inner) {
			return false;
		}
		/// points order does matter
		for (int i = 0; i < N; i++) {
			if (this->p[i] != other.p[i]) {
				return false;
			}
		}
		return true;
	}
	
	bool operator!=(const Element& other) const {
		return !( (*this) == other );
	}
	
	/** Return all points owned by both elements */
	std::set<Point> equalPoints(const Element& other) {
		assert_true(this->inner && other.inner);
		std::set<Point> ans;
		
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (this->p[i] == other.p[j]) {
					ans.insert(this->p[i]);
					break;
				}
			}
		}
		
		return ans;
	}
	
};

template<typename Point> using Triangle = Element<Point, 3>;

}

}

#endif // LIBGCM_ELEMENTS_HPP
