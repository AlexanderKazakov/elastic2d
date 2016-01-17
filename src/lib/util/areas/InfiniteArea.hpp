#ifndef LIBGCM_INFINITEAREA_HPP
#define LIBGCM_INFINITEAREA_HPP

#include "lib/util/areas/Area.hpp"

namespace gcm {
	/**
	 * Area that contains any Vector
	 */
	class InfiniteArea : public Area {
	public:
		InfiniteArea() { };
		virtual bool contains(const linal::Vector3& coords) const override {
			return true;
		};
	};
}

#endif // LIBGCM_INFINITEAREA_HPP
