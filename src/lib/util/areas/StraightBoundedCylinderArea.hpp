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
		StraightBoundedCylinderArea(const real& _radius, const linal::Vector3& _begin, const linal::Vector3& _end);
		virtual bool contains(const linal::Vector3& coords) const override;
		virtual void move(const linal::Vector3& shift) override;
	private:
		real radius;
		linal::Vector3 begin; // center of upper cap of the keg
		linal::Vector3 end; // center of lower cap of the keg
		linal::Vector3 axis; // normalized vector along axis of the cylinder
	};
}

#endif // LIBGCM_CYLINDERAREA_HPP
