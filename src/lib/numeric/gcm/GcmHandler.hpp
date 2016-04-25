#ifndef LIBGCM_GCMHANDLER_HPP
#define LIBGCM_GCMHANDLER_HPP

#include <lib/linal/linal.hpp>
#include <lib/mesh/DefaultMesh.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>
#include <lib/mesh/grid/Cgal2DGrid.hpp>
#include <lib/numeric/interpolation/interpolation.hpp>


namespace gcm {
/**
 * Bridge between grids of different types and grid-characteristic method
 */
template<typename TModel, typename TGrid, typename TMaterial>
struct GcmHandler;

template<typename TModel, typename TMaterial>
struct GcmHandler<TModel, CubicGrid, TMaterial> {
	typedef CubicGrid                            Grid;
	typedef DefaultMesh<TModel, Grid, TMaterial> Mesh;
	typedef typename Mesh::Matrix                Matrix;
	typedef typename Mesh::PdeVector             PdeVector;
	typedef typename Mesh::Iterator              Iterator;

	/**
	 * Interpolate nodal values in specified points.
	 * Interpolated value for k-th point in vector dx
	 * stored in k-th column of returned Matrix.
	 * @param stage direction
	 * @param it index-iterator of node
	 * @param dx Vector of distances from reference node on which
	 * values should be interpolated
	 * @return Matrix with interpolated nodal values in columns
	 */
	Matrix interpolateValuesAround(const Mesh& mesh, const int stage,
	                               const Iterator& it, const PdeVector& dx) const {
		Matrix ans;
		std::vector<PdeVector> src( (size_t) (mesh.borderSize + 1) );
		PdeVector res;
		Iterator shift({0, 0, 0}, mesh.sizes);
		for (int k = 0; k < PdeVector::M; k++) {
			shift(stage) = (dx(k) > 0) ? 1 : -1;
			for (int i = 0; i < src.size(); i++) {
				src[(size_t)i] = mesh.pde(it + shift * i);
			}
			EqualDistanceLineInterpolator<PdeVector>::minMaxInterpolate(
					res, src, fabs(dx(k)) / mesh.h(stage));
			ans.setColumn(k, res);
		}
		return ans;
	}

	void borderCorrector(const Mesh&, const int, const Iterator&) const { }
	void beforeStage(const Mesh&) { }
};

template<typename TModel, typename TMaterial>
struct GcmHandler<TModel, Cgal2DGrid, TMaterial> {
	typedef Cgal2DGrid                              Grid;
	typedef TModel                                  Model;
	typedef DefaultMesh<Model, Grid, TMaterial>     Mesh;
	typedef typename Mesh::Matrix                   Matrix;
	typedef typename Mesh::PdeVector                PdeVector;
	typedef typename Mesh::Iterator                 Iterator;
	typedef typename Mesh::BORDER_CONDITION         BORDER_CONDITION;
	typedef linal::VECTOR<2, PdeVector>             PdeGradient;
	typedef linal::SYMMETRIC_MATRIX<2, PdeVector>   PdeHessian;
	
	static const int OUTER_NUMBER = BORDER_CONDITION::OUTER_NUMBER;
	static const int PDE_SIZE = Model::PDE_SIZE;
	
	typedef linal::Matrix<PDE_SIZE, OUTER_NUMBER> OuterU1Matrix;

	/**
	 * Interpolate nodal values in specified points.
	 * Interpolated value for k-th point in vector dx are
	 * stored in k-th column of returned Matrix.
	 * If specified point appears out of body, matrix column is set to zeros
	 * and outerInvariants is added with the index.
	 * @param stage direction
	 * @param it index-iterator of node
	 * @param dx Vector of distances from reference node on which
	 * values should be interpolated
	 * @return Matrix with interpolated nodal values in columns
	 */
	Matrix interpolateValuesAround(const Mesh& mesh, const int stage,
	                               const Iterator& it, const PdeVector& dx) {
	    outerInvariants.clear();
		Matrix ans;
		
		for (int k = 0; k < PdeVector::M; k++) {
			
			if (dx(k) == 0) {
			// special for exact hit
				ans.setColumn(k, mesh.pde(it));
				continue;
			}
			
			// point to interpolate respectively to point by iterator it
			Real2 shift({(stage == 0) * dx(k), (stage == 1) * dx(k)});
						
			auto t = mesh.findOwnerTriangle(it, shift);
			PdeVector u; linal::clear(u);
			if (t.inner) {
			// interpolate inner value
				u = TriangleInterpolator<PdeVector>::interpolate(
					mesh.coords2d(t.p[0]), mesh.pde(t.p[0]), gradients[t.p[0].iter],
					mesh.coords2d(t.p[1]), mesh.pde(t.p[1]), gradients[t.p[1].iter],
					mesh.coords2d(t.p[2]), mesh.pde(t.p[2]), gradients[t.p[2].iter],
					mesh.coords2d(it) + shift);
			} else {
			// add outer invariant for border corrector
				outerInvariants.push_back(k);
			}
			ans.setColumn(k, u);
		}
		
		return ans;
	}
	
	
	/**
	 * Apply border conditions according to Chelnokov's PhD thesis, page 42.
	 * Here used outerInvariants, written at interpolateValuesAround before.
	 */
	void borderCorrector(Mesh& mesh, const int stage, const Iterator& it) {
		if (outerInvariants.size() == 0) { return; }
		
		if (mesh.borderIndices.find(it.iter) == mesh.borderIndices.end()) {
			LOG_DEBUG("Bad case: " << outerInvariants.size() << " outer invariants "
				<< "in inner node at " << mesh.coords2d(it));
			return;
		}
		
		if (outerInvariants.size() != OUTER_NUMBER) {
			LOG_DEBUG("Bad case: " << outerInvariants.size() << " outer invariants "
				<< "in border node at " << mesh.coords2d(it));
			return;
		}
		
		const BORDER_CONDITION* borderCondition = mesh.getBorderCondition(it);
		if (borderCondition == nullptr) { return; } // non-reflection
		
		const Real2 normal = mesh.normal(it);
		auto B = borderCondition->B(normal);
		auto b = borderCondition->b(normal);
		OuterU1Matrix outerU1;	
		for (int i = 0; i < OUTER_NUMBER; i++) {
			outerU1.setColumn(
				i, mesh.matrices(it)->m[stage].U1.getColumn(outerInvariants[(size_t)i]));
		}
		
		const auto alpha = linal::solveLinearSystem(B * outerU1, b - B * mesh.pdeNew(it));
		mesh._pdeNew(it) = mesh.pdeNew(it) + outerU1 * alpha;
	}
	
	
	void beforeStage(const Mesh& mesh) {
		/// calculate spatial derivatives of all mesh pde values ones before stage
		/// in order to use them multiple times while stage calculation 
		estimateGradient(mesh);
	}
	

private:
	/// List of outer Riemann invariants used in borderCorrector.
	/// Invariants are specified by their indices in matrix L.
	std::vector<int> outerInvariants;
	
	/// The storage of gradients of mesh pde values.
	std::vector<PdeGradient> gradients;
	/// The storage of hessians of mesh pde values. (Unused now)
	std::vector<PdeHessian> hessians;
	
	USE_AND_INIT_LOGGER("gcm.Cgal2DGridGcmHandler")
	
	
	void estimateGradient(const Mesh& mesh) {
		gradients.resize(mesh.sizeOfAllNodes());

		for (auto v : mesh) {
			/// number of neighbor vertices in 2D is unlikely to be more than 8			
			auto A = linal::Matrix<8, 2>::Zeros(); // SLE matrix
			auto b = linal::VECTOR<8, PdeVector>::Zeros(); // SLE right part
			auto W = linal::DiagonalMatrix<8>::Zeros(); // weight matrix for
					// linear least squares method
			
			/// estimate gradient from its definition:
			/// \f$  (\nabla f, \vec{a} - \vec{b}) = f(a) - f(b),  \f$ 
			/// applying it to each neighbor of the vertex
			auto neighbors = mesh.findNeighborVertices(v);
			int i = 0;
			for (auto neighbor : neighbors) {
				Real2 d = mesh.coords2d(neighbor) - mesh.coords2d(v);
				A.setRow(i, d);
				W(i) = linal::length(d);
				b(i) = mesh.pde(neighbor) - mesh.pde(v);
				
				i++; if (i > 7) { break; }
			}
			
			gradients[v.iter] = linal::linearLeastSquares(A, b, W);
		}
	}
	
	
	void estimateHessian(const Mesh& mesh) {
	/// the function is currently unused because third order is not 
	/// so straightforward as second order
		hessians.resize(mesh.sizeOfAllNodes());
		assert_eq(hessians.size(), gradients.size());

		for (auto v : mesh) {
			/// number of neighbor vertices in 2D is unlikely to be more than 8			
			auto A = linal::Matrix<8, 2>::Zeros();
			auto b = linal::VECTOR<8, PdeGradient>::Zeros();
			
			/// estimate Hessian as the gradient of gradient:
			/// \f$  \matrix{H}(f) = \nabla (\nabla f),  \f$ 
			/// applying it to each neighbor of the vertex
			auto neighbors = mesh.findNeighborVertices(v);
			int i = 0;
			for (auto neighbor : neighbors) {
				A.setRow(i, mesh.coords2d(neighbor) - mesh.coords2d(v));
				b(i) = gradients[neighbor.iter] - gradients[v.iter];
				
				i++; if (i > 7) { break; }
			}
			
			auto H = linal::linearLeastSquares(A, b); // gradient of gradient
			
			hessians[v.iter] = {H(0)(0), (H(0)(1) + H(1)(0)) / 2.0, H(1)(1)};
		}
	}
};


}

#endif // LIBGCM_GCMHANDLER_HPP
