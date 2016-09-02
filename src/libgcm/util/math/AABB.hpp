#ifndef LIBGCM_AXESALIGNEDBOUNDARYBOX_HPP
#define LIBGCM_AXESALIGNEDBOUNDARYBOX_HPP

#include <libgcm/util/infrastructure/infrastructure.hpp>

namespace gcm {


/**
 * Axes-aligned boundary box (AABB)
 */
template<typename Point>
struct AxesAlignedBoundaryBox {
	static const int D = Point::SIZE;
	
	/// the most left point of the box
	Point min;
	/// the most right point of the box
	Point max;
	
	
	/**
	 * AABB includes its min and max points,
	 * so the case (min == max) is valid
	 */
	bool valid() const {
		for (int i = 0; i < D; i++) {
			if (sizes()(i) < 0) { return false; }
		}
		return true;
	}
	
	
	/**
	 * The AABB is "slice" if its length along some direction is zero.
	 */
	int sliceDirection() const {
		for (int i = 0; i < D; i++) {
			if (sizes()(i) == 0) { return i; }
		}
		THROW_INVALID_ARG("The AABB has non-zero lengthes along all directions");
	}
	
	
	/**
	 * Lengthes of the AABB in each direction
	 */
	Point sizes() const {
		return max - min;
	}
	
	
	/**
	 * Move AABB on specified direction
	 */
	static AxesAlignedBoundaryBox translate(
			const AxesAlignedBoundaryBox& b, const Point& direction) {
		return {b.min + direction, b.max + direction};
	}
	
	
	/**
	 * Intersection of two AABBs.
	 * If there is no intersection along some axis,
	 * intersection.min will be greater than intersection.max,
	 * i.e. invalid.
	 */
	static AxesAlignedBoundaryBox intersection(
			const AxesAlignedBoundaryBox& a, const AxesAlignedBoundaryBox& b) {
		AxesAlignedBoundaryBox ans;
		for (int i = 0; i < D; i++) {
			ans.min(i) = std::max(a.min(i), b.min(i));
			ans.max(i) = std::min(a.max(i), b.max(i));
		}
		return ans;
	}
	
};


}


namespace std {

template<typename Point>
inline std::ostream& operator<<(std::ostream& os,
		const gcm::AxesAlignedBoundaryBox<Point>& aabb) {
	
	os << "\nmin: " << aabb.min << "\nmax: " << aabb.max;
	
	return os;
}

}


#endif // LIBGCM_AXESALIGNEDBOUNDARYBOX_HPP
