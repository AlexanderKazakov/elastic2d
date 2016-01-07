#ifndef LIBGCM_IDEALELASTIC3DMODEL_HPP
#define LIBGCM_IDEALELASTIC3DMODEL_HPP

#include "lib/model/Model.hpp"
#include "lib/nodes/IdealElastic3DNode.hpp"
#include "lib/gcm_matrices/IdealElastic3DGcmMatrices.hpp"

namespace gcm {
	class IdealElastic3DModel : public Model {
	public:
		typedef IdealElastic3DNode Node;
		typedef IdealElastic3DGcmMatrices GcmMatrices;
		static const int DIMENSIONALITY = GcmMatrices::DIMENSIONALITY;
	};
}


#endif // LIBGCM_IDEALELASTIC3DMODEL_HPP
