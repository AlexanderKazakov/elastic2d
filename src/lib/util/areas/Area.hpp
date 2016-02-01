#ifndef LIBGCM_AREA_HPP
#define LIBGCM_AREA_HPP

#include <lib/linal/special/Vector3.hpp>
#include <lib/linal/LinalRoutines.hpp>

namespace gcm {
	class Area {
	public:
		/**
		 * @return true if the area contains specified point, false otherwise
		 * (by definition, area does not contain points on its border)
		 */
		virtual bool contains(const linal::Vector3& coords) const = 0;
	};
}

#endif // LIBGCM_AREA_HPP
