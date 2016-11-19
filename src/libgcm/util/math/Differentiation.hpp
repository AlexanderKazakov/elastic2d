#ifndef LIBGCM_DIFFERENTIATION_HPP
#define LIBGCM_DIFFERENTIATION_HPP

#include <libgcm/linal/linal.hpp>

namespace gcm {

/**
 * Class for spatial differentiation. Given with some mesh (function data
 * and geometry), calculate gradient, hessian, etc.
 */
template<typename TMesh>
class Differentiation {
public:
	typedef TMesh                                              Mesh;
	typedef typename Mesh::Grid                                Grid;
	typedef typename Mesh::PdeVector                           PdeVector;

	static const int DIMENSIONALITY = Grid::DIMENSIONALITY;
	
	typedef linal::VECTOR<DIMENSIONALITY, PdeVector>           PdeGradient;
	typedef linal::SYMMETRIC_MATRIX<DIMENSIONALITY, PdeVector> PdeHessian;
	typedef linal::Vector<DIMENSIONALITY>                      RealD;
	
	static const int MAX_NUMBER_OF_NEIGHBOR_VERTICES = 
			Grid::MAX_NUMBER_OF_NEIGHBOR_VERTICES;
	
	
	/** 
	 * Calculate gradients of mesh pde values at each mesh vertex.
	 * The order in answer is equal to order of values in the mesh.
	 */
	static void estimateGradient(const Mesh& mesh,
			std::vector<PdeGradient>& gradients) {
			
		gradients.resize(mesh.sizeOfAllNodes());
		
//		#pragma omp parallel for
		for (size_t it = 0; it < mesh.sizeOfRealNodes(); ++it) {
			// SLE matrix
			auto A = linal::Matrix<MAX_NUMBER_OF_NEIGHBOR_VERTICES, DIMENSIONALITY>::Zeros();
			// SLE right part
			auto b = linal::VECTOR<MAX_NUMBER_OF_NEIGHBOR_VERTICES, PdeVector>::Zeros();
			// weight matrix for linear least squares method
			auto W = linal::DiagonalMatrix<MAX_NUMBER_OF_NEIGHBOR_VERTICES>::Zeros();
			
			/// estimate gradient from its definition:
			/// \f$  (\nabla f, \vec{a} - \vec{b}) = f(a) - f(b),  \f$ 
			/// applying it to each neighbor of the vertex
			const auto neighbors = mesh.findNeighborVertices(it);
			int i = 0;
			for (const auto neighbor : neighbors) {
				RealD d = mesh.coordsD(neighbor) - mesh.coordsD(it);
				A.setRow(i, d);
				W(i) = 1.0 / linal::length(d);
				b(i) = mesh.pde(neighbor) - mesh.pde(it);
				
				i++; if (i == MAX_NUMBER_OF_NEIGHBOR_VERTICES) { break; }
			}
			
			gradients[mesh.getIndex(it)] = linal::linearLeastSquares(A, b, W);
		}
	}
	
	
	/** 
	 * Calculate Hessians of mesh pde values at each mesh vertex.
	 * The order of gradients must be equal to order of values in mesh.
	 * The order in answer is equal to order of values in the mesh.
	 */
	static void estimateHessian(const Mesh& mesh, 
			const std::vector<PdeGradient>& gradients,
			std::vector<PdeHessian>& hessians) {
		
		hessians.resize(mesh.sizeOfAllNodes());
		assert_eq(hessians.size(), gradients.size());

		// SLE matrix
		linal::Matrix<MAX_NUMBER_OF_NEIGHBOR_VERTICES, DIMENSIONALITY> A;
		// SLE right part
		linal::VECTOR<MAX_NUMBER_OF_NEIGHBOR_VERTICES, PdeGradient> b;
		// weight matrix for linear least squares method
		linal::DiagonalMatrix<MAX_NUMBER_OF_NEIGHBOR_VERTICES> W;

		for (const auto v : mesh) {
			linal::clear(A); linal::clear(b); linal::clear(W);
			
			/// estimate Hessian as the gradient of gradient:
			/// \f$  \matrix{H}(f) = \nabla (\nabla f),  \f$ 
			/// applying it to each neighbor of the vertex
			const auto neighbors = mesh.findNeighborVertices(v);
			int i = 0;
			for (const auto neighbor : neighbors) {
				RealD d = mesh.coordsD(neighbor) - mesh.coordsD(v);
				A.setRow(i, d);
				W(i) = 1.0 / linal::length(d);
				b(i) = gradients[mesh.getIndex(neighbor)] - gradients[mesh.getIndex(v)];
				
				i++; if (i == MAX_NUMBER_OF_NEIGHBOR_VERTICES) { break; }
			}
			
			auto H = linal::linearLeastSquares(A, b, W); // gradient of gradient
			
			hessians[mesh.getIndex(v)] = {H(0)(0), (H(0)(1) + H(1)(0)) / 2.0, H(1)(1)};
		}
	}
	
};


}

#endif // LIBGCM_DIFFERENTIATION_HPP
