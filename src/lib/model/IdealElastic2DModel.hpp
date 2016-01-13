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

		static const std::map<const std::string /* name */,
				const std::pair<const int /* index */, const int /* size */>> VECTORS;

		static const std::map<const std::string /* name */,
				const int /* index */> SCALARS;

	};
}


#endif // LIBGCM_IDEALELASTIC2DMODEL_HPP
