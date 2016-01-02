#ifndef LIBGCM_IDEALELASTIC2DMODEL_HPP
#define LIBGCM_IDEALELASTIC2DMODEL_HPP

#include "lib/model/Model.hpp"
#include "lib/nodes/IdealElastic2DNode.hpp"

namespace gcm {
	class IdealElastic2DModel : public Model {
		typedef IdealElastic2DNode Node;

	};
}


#endif // LIBGCM_IDEALELASTIC2DMODEL_HPP
