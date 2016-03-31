#ifndef LIBGCM_CYLINDERAREA_HPP
#define LIBGCM_CYLINDERAREA_HPP

#include <lib/util/areas/Area.hpp>

namespace gcm {
/** Straight bounded cylinder (aka barrel, keg, etc..) */
class StraightBoundedCylinderArea : public Area {
public:
	/**
	 * Straight bounded cylinder (aka barrel, keg, etc..)
	 * @param radius radius
	 * @param begin center of upper cap of the keg
	 * @param end center of lower cap of the keg
	 */
	StraightBoundedCylinderArea(const real& _radius, const Real3& _begin, const Real3& _end);
	virtual bool contains(const Real3& coords) const override;
	virtual void move(const Real3& shift) override;

private:
	real radius;
	Real3 begin;       ///< center of upper cap of the keg
	Real3 end;         ///< center of lower cap of the keg
	Real3 axis;        ///< normalized vector along axis of the cylinder
};


}

#endif // LIBGCM_CYLINDERAREA_HPP
