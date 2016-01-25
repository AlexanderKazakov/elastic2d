#ifndef LIBGCM_SPHEREAREA_HPP
#define LIBGCM_SPHEREAREA_HPP

#include <lib/util/areas/Area.hpp>

namespace gcm {
	class SphereArea : public Area {
	public:
		SphereArea(const real& _radius, const linal::Vector3& _center);
		virtual bool contains(const linal::Vector3& coords) const override;
	private:
		real radius;
		linal::Vector3 center;
	};
}

#endif // LIBGCM_SPHEREAREA_HPP
