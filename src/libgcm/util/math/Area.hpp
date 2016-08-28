#ifndef LIBGCM_AREA_HPP
#define LIBGCM_AREA_HPP

#include <libgcm/linal/linal.hpp>

namespace gcm {
struct Area {
	/**
	 * @return true if the area contains specified point, false otherwise
	 * (by definition, area does not contain points on its border)
	 */
	virtual bool contains(const Real3& coords) const = 0;

	/**
	 * Move area on specified distances
	 */
	virtual void move(const Real3& shift) = 0;

};


struct InfiniteArea : public Area {
	virtual bool contains(const Real3&) const override { return true; }
	virtual void move(const Real3&) { }
};


struct AxisAlignedBoxArea : public Area {
	/**
	 * @param min the nearest point to origin of coordinates
	 * @param max the farthest point from origin of coordinates
	 */
	AxisAlignedBoxArea(const Real3& _min, const Real3& _max) :
			min(_min), max(_max) {
		for (int i = 0; i < 3; i++) {
			assert_gt((max - min)(i), 0.0);
		}
	}
	
	virtual bool contains(const Real3& coords) const override {
		for (int i = 0; i < 3; i++) {
			if (coords(i) <= min(i) || coords(i) >= max(i)) {return false; }
		}
		return true;
	}
	
	virtual void move(const Real3& shift) override {
		min += shift; max += shift;
	}

	Real3 getMin() const { return min; }
	Real3 getMax() const { return max; }

private:
	Real3 min;         ///< the nearest point to origin of coordinates
	Real3 max;         ///< the farthest point from origin of coordinates
};


struct SphereArea : public Area {
	SphereArea(const real& _radius, const Real3& _center) :
			radius(_radius), center(_center) {
		assert_gt(radius, 0.0);
	}
	
	virtual bool contains(const Real3& coords) const override {
		return linal::length(coords - center) < radius;
	}
	
	virtual void move(const Real3& shift) override {
		center += shift;
	}
	
private:
	real radius;
	Real3 center;
};


/** Straight bounded cylinder (aka barrel, keg, etc..) */
class StraightBoundedCylinderArea : public Area {
public:
	/**
	 * Straight bounded cylinder (aka barrel, keg, etc..)
	 * @param radius radius
	 * @param begin center of upper cap of the keg
	 * @param end center of lower cap of the keg
	 */
	StraightBoundedCylinderArea(const real& _radius,
			const Real3& _begin, const Real3& _end) :
		radius(_radius), begin(_begin), end(_end) {
		
		axis = linal::normalize(end - begin);
		assert_gt(radius, 0.0);
	}
	
	virtual bool contains(const Real3& coords) const override {
		/// check that the point is between caps
		if (linal::dotProduct(coords - begin, axis) *
			linal::dotProduct(coords - end, axis) >= 0) { 
			return false;
		}
	
		/// check that the point is on the less than radius distance from the axis
		real projection = linal::dotProduct(coords - begin, axis);
		return linal::dotProduct(coords - begin, coords - begin) -
			projection * projection < radius * radius;
	}
	
	virtual void move(const Real3& shift) override {
		begin += shift; end += shift;
	}

private:
	real radius;       ///< radius of cylinder
	Real3 begin;       ///< center of upper cap of the keg
	Real3 end;         ///< center of lower cap of the keg
	Real3 axis;        ///< normalized vector along axis of the cylinder
};


}

#endif // LIBGCM_AREA_HPP
