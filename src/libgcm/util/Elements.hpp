#ifndef LIBGCM_ELEMENTS_HPP
#define LIBGCM_ELEMENTS_HPP

#include <set>

#include <libgcm/linal/linal.hpp>


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
	int n = 0;               ///< number of set points, 0 <= n <= N
	
	Element() = default;
	
	template<typename ForwIter, typename Predicate>
	Element(ForwIter begin, const ForwIter end, const Predicate transformFunc) {
		n = 0;
		for (; begin != end; ++begin) {
			assert_lt(n, N);
			(*this)(n++) = transformFunc(*begin);
		}
	}
	
	Element(std::initializer_list<Point> list) :
			Element(list.begin(), list.end()) { }
	
	template<typename OtherPointType, typename Predicate>
	Element(const Element<OtherPointType, N>& other,
			const Predicate transformFunc) {
		this->n = other.n;
		for (int i = 0; i < n; i++) {
			this->p[i] = transformFunc(other(i));
		}
	}
	
	const Point& operator()(const int i) const { return p[i]; }
	      Point& operator()(const int i)       { return p[i]; }
	
	bool operator==(const Element& other) const {
		if (this->n != other.n) {
			return false;
		}
		/// points order does matter
		for (int i = 0; i < n; i++) {
			if (this->p[i] != other(i)) {
				return false;
			}
		}
		return true;
	}
	
	bool operator!=(const Element& other) const {
		return !( (*this) == other );
	}
	
	/** 
	 * Return all points owned by both elements.
	 * Elements must have equal number of set points.
	 */
	std::set<Point> equalPoints(const Element& other) const {
		assert_eq(this->n, other.n);
		std::set<Point> ans;
		
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (this->p[i] == other(j)) {
					ans.insert(this->p[i]);
					break;
				}
			}
		}
		
		return ans;
	}
	
	
	bool has(const Point& point) const {
		for (int i = 0; i < n; i++) {
			if ((*this)(i) == point) {
				return true;
			}
		}
		return false;
	}
	
};


template<typename Point, int TN>
std::vector<Point> sortedUniquePointsOfElements(
		const std::vector<Element<Point, TN>>& items) {
	std::set<Point> points;
	for (const Element<Point, TN>& item : items) {
		assert_eq(item.n, item.N);
		for (int i = 0; i < item.n; i++) {
			points.insert(item(i));
		}
	}
	return std::vector<Point>(points.begin(), points.end());
}


template<typename Point> using Triangle = Element<Point, 3>;
template<typename Point> using Tetrahedron = Element<Point, 4>;


}

}


namespace std {

template<typename Point, int TN>
inline std::ostream& operator<<(std::ostream& os, 
		const gcm::elements::Element<Point, TN>& element) {

	os << "Number of set points (n) : " << element.n << std::endl;
	if (element.n > 0) {
		os << "Points: " << std::endl;
		for (int j = 0; j < element.n; j++) {
			os << element(j) << "\t";
		}
		os << std::endl;
	}
	return os;
}


}

#endif // LIBGCM_ELEMENTS_HPP
