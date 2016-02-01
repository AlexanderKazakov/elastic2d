#ifndef LIBGCM_NODE_HPP
#define LIBGCM_NODE_HPP

#include <memory>

namespace gcm {
	/* ====================================================================

	  This is building blocks for grids :
	  - simple nodes for homogeneous linear models and structured grids
	  - nodes with pointer to rheology matrices for inhomogeneous or nonlinear
	  models and structured grids
	  - nodes with coordinates and numbers for unstructured grids
	  - cells and surfaces
	  - etc...

	  The blocks are constructed by multiple inheritance from special
	  "Has/HasNot" storage classes.
	  The reason for inheritance not composition see at lib/test/sequence/Demo.cpp

	==================================================================== */

	template<typename TModel>
	struct HasGcmMatrices {
		std::shared_ptr<typename TModel::GCM_MATRICES> matrix;
	};
	template<typename TModel> struct HasNotGcmMatrices { };
	template<typename TModel> using DefaultGcmMatricesStorage = HasGcmMatrices<TModel>;

	template<typename TModel,
			template<typename> class GcmMatricesStorage = DefaultGcmMatricesStorage>
	struct Node : public TModel::OdeVariables, public GcmMatricesStorage<TModel> {

		typedef TModel Model;
		typedef typename Model::Variables Variables;
		typedef typename Model::Vector Vector;
		typedef typename Model::GCM_MATRICES GCM_MATRICES;

		static MPI::Datatype MPI_NODE_TYPE; // Special type for node for MPI connection

		Vector u; // vector of PDE variables
	};

	template<typename TModel, template<typename> class GcmMatricesStorage>
	MPI::Datatype Node<TModel, GcmMatricesStorage>::MPI_NODE_TYPE;
}

#endif // LIBGCM_NODE_HPP
