#ifndef LIBGCM_IDEALELASTIC2DMODEL_HPP
#define LIBGCM_IDEALELASTIC2DMODEL_HPP

#include "lib/model/Model.hpp"
#include "lib/nodes/IdealElastic2DNode.hpp"
#include "lib/gcm_matrices/IdealElastic2DGcmMatrices.hpp"

namespace gcm {
	class IdealElastic2DModel : public Model {
	public:
		typedef IdealElastic2DNode Node;
		typedef IdealElastic2DGcmMatrices GcmMatrices;
		static const int DIMENSIONALITY = GcmMatrices::DIMENSIONALITY;
	};
}


#endif // LIBGCM_IDEALELASTIC2DMODEL_HPP
