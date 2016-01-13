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

		static const std::map<const std::string /* name */,
				const std::pair<const int /* index */, const int /* size */>> VECTORS;

		static const std::map<const std::string /* name */,
				const int /* index */> SCALARS;

	};
}


#endif // LIBGCM_IDEALELASTIC3DMODEL_HPP
