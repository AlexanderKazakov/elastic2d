#ifndef LIBGCM_IDEALELASTIC1DMODEL_HPP
#define LIBGCM_IDEALELASTIC1DMODEL_HPP

#include "lib/model/Model.hpp"
#include "lib/nodes/IdealElastic1DNode.hpp"
#include "lib/gcm_matrices/IdealElastic1DGcmMatrices.hpp"

namespace gcm {
	class IdealElastic1DModel : public Model {
	public:
		typedef IdealElastic1DNode Node;
		typedef IdealElastic1DGcmMatrices GcmMatrices;
		static const int DIMENSIONALITY = GcmMatrices::DIMENSIONALITY;
	};
}


#endif // LIBGCM_IDEALELASTIC1DMODEL_HPP
