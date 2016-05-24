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

	static const int N = TN; ///< number of points
	Point p[N];              ///< points
	bool valid = false;      ///< indicates that points are set
	
	Element() = default;
	
	Element(std::initializer_list<Point> list) {
		assert_eq(list.size(), N);
		std::copy(list.begin(), list.end(), p);
	}
	
	template<typename OtherPointType>
	Element(const Element<OtherPointType, N>& other,
			std::function<Point(const OtherPointType&)> transformFunc) {
		this->valid = other.valid;
		for (int i = 0; i < N; i++) {
			this->p[i] = transformFunc(other(i));
		}
	}
	
	const Point& operator()(const int i) const { return p[i]; }
	      Point& operator()(const int i)       { return p[i]; }
	
	bool operator==(const Element& other) const {
		assert_true(this->valid && other.valid);
		/// points order does matter
		for (int i = 0; i < N; i++) {
			if (this->p[i] != other(i)) {
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
		assert_true(this->valid && other.valid);
		std::set<Point> ans;
		
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (this->p[i] == other(j)) {
					ans.insert(this->p[i]);
					break;
				}
			}
		}
		
		return ans;
	}
	
};


template<typename Point> using Triangle = Element<Point, 3>;
template<typename Point> using Tetrahedron = Element<Point, 4>;


}

}


namespace std {

template<typename Point, int TN>
inline std::ostream& operator<<(std::ostream& os, 
		const gcm::elements::Element<Point, TN>& element) {

	os << "Element valid : " << element.valid << std::endl;
	if (element.valid) {
		os << "Points: " << std::endl;
		for (int j = 0; j < TN; j++) {
			os << element(j) << "\t";
		}
		os << std::endl;
	}
	return os;
}


}

#endif // LIBGCM_ELEMENTS_HPP
