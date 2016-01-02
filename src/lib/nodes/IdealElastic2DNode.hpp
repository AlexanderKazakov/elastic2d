#ifndef LIBGCM_IDEALELASTIC2DNODE_HPP
#define LIBGCM_IDEALELASTIC2DNODE_HPP

#include <string.h>
#include "lib/nodes/Node.hpp"

namespace gcm {
	class IdealElastic2DNode : public Node<5, IdealElastic2DContainer> {

	};

	class IdealElastic2DContainer {
	public:
		static const int SIZE = 5;
		union {
			real values[SIZE];
			struct {
				real Vx;
				real Vy;
				real Sxx;
				real Sxy;
				real Syy;
			};
		};
	};
}


#endif // LIBGCM_IDEALELASTIC2DNODE_HPP
