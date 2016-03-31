#ifndef LIBGCM_SPHEREAREA_HPP
#define LIBGCM_SPHEREAREA_HPP

#include <lib/util/areas/Area.hpp>

namespace gcm {
class SphereArea : public Area {
public:
	SphereArea(const real& _radius, const Real3& _center);
	virtual bool contains(const Real3& coords) const override;
	virtual void move(const Real3& shift) override;

private:
	real radius;
	Real3 center;
};


}

#endif // LIBGCM_SPHEREAREA_HPP
