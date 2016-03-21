#ifndef LIBGCM_BOXAREA_HPP
#define LIBGCM_BOXAREA_HPP

#include <lib/util/areas/Area.hpp>

namespace gcm {
	class AxisAlignedBoxArea : public Area {
	public:
		/**
		 * @param min the nearest point to origin of coordinates
		 * @param max the farthest point from origin of coordinates
		 */
		AxisAlignedBoxArea(const Real3& _min, const Real3& _max);
		virtual bool contains(const Real3& coords) const override;
		virtual void move(const Real3& shift) override;
		Real3 getMin() const { return min; };
		Real3 getMax() const { return max; };
	private:
		Real3 min; // the nearest point to origin of coordinates
		Real3 max; // the farthest point from origin of coordinates
	};
}

#endif // LIBGCM_BOXAREA_HPP
