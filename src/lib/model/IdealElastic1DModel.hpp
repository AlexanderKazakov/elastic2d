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

		static const std::map<const std::string /* name */,
				const std::pair<const int /* index */, const int /* size */>> VECTORS;

		static const std::map<const std::string /* name */,
				const int /* index */> SCALARS;

	};
}


#endif // LIBGCM_IDEALELASTIC1DMODEL_HPP
