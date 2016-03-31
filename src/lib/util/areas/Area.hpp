#ifndef LIBGCM_AREA_HPP
#define LIBGCM_AREA_HPP

#include <lib/linal/linal.hpp>

namespace gcm {
class Area {
public:
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


}

#endif // LIBGCM_AREA_HPP
