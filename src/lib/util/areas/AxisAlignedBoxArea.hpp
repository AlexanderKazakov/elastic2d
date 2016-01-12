#ifndef LIBGCM_BOXAREA_HPP
#define LIBGCM_BOXAREA_HPP

#include "lib/util/areas/Area.hpp"

namespace gcm {
	class AxisAlignedBoxArea : public Area {
	public:
		/**
		 * @param min the nearest point to origin of coordinates
		 * @param max the farthest point from origin of coordinates
		 */
		AxisAlignedBoxArea(const linal::Vector3& min, const linal::Vector3& max);
		virtual bool contains(const linal::Vector3& coords) const override;
	private:
		linal::Vector3 min; // the nearest point to origin of coordinates
		linal::Vector3 max; // the farthest point from origin of coordinates
	};
}

#endif // LIBGCM_BOXAREA_HPP
