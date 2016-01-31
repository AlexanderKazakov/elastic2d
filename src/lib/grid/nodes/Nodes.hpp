#ifndef LIBGCM_NODE_HPP
#define LIBGCM_NODE_HPP

#include <memory>

namespace gcm {
	/* ==================================================================
	  This is building blocks for grids :
	  - simple nodes for homogeneous linear models and structured grids
	  - nodes with pointer to rheology matrices for inhomogeneous or nonlinear
	  models and structured grids
	  - nodes with coordinates and numbers for unstructured grids
	  - cells and surfaces
	  - etc..
	   ================================================================== */

	/**
	 * Simple node contains only vector of PDE variables
	 */
	template<typename TModel>
	struct Node {
		typedef typename TModel::Variables Variables;
		typedef typename TModel::Vector Vector;
		typedef typename TModel::GCM_MATRICES GCM_MATRICES;

		// TODO - for that type we can use byte copy, isn't it?
		// static MPI::Datatype MPI_NODE_TYPE; // Special type for node for MPI connection

		Vector u; // vector of PDE variables
	};

	/**
	 * Node with vector of PDE variables and pointer to rheology matrices
	 */
	template<typename TModel>
	struct NodeMatrix {
		typedef typename TModel::Variables Variables;
		typedef typename TModel::Vector Vector;
		typedef typename TModel::GCM_MATRICES GCM_MATRICES;

		static MPI::Datatype MPI_NODE_TYPE; // Special type for node for MPI connection

		std::shared_ptr<GCM_MATRICES> matrix;
		Vector u; // vector of PDE variables
	};

	template<class TModel> MPI::Datatype NodeMatrix<TModel>::MPI_NODE_TYPE;
}

#endif // LIBGCM_NODE_HPP
