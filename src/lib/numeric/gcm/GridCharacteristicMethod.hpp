#ifndef LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
#define LIBGCM_GRIDCHARACTERISTICMETHOD_HPP

#include <lib/linal/linal.hpp>
#include <lib/mesh/DefaultMesh.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>
#include <lib/mesh/grid/Cgal2DGrid.hpp>
#include <lib/numeric/interpolation/interpolation.hpp>


namespace gcm {

/** @return values on next time layer by grid-characteristic method */
template<typename TMatrix>
inline auto
localGcmStep(const TMatrix& U1, const TMatrix& U, const TMatrix& values) 
		-> decltype(U1 * linal::diagonalMultiply(U, values)) {
	return /* new values = U1 * Riemann solvers */ U1 * linal::diagonalMultiply(
	       /* Riemann solvers = U * old values */ U,
	       /* old values are in columns of the matrix */ values);
}


/**
 * Grid-characteristic method
 */
template<typename TModel, typename TGrid, typename TMaterial>
class GridCharacteristicMethod;


/**
 * Grid-characteristic method specialization for cubic grids
 */
template<typename TModel, typename TMaterial>
class GridCharacteristicMethod<TModel, CubicGrid, TMaterial> {
public:
	typedef CubicGrid                            Grid;
	typedef DefaultMesh<TModel, Grid, TMaterial> Mesh;
	typedef typename Mesh::Matrix                Matrix;
	typedef typename Mesh::PdeVector             PdeVector;
	typedef typename Mesh::Iterator              Iterator;

	/**
	 * Do grid-characteristic stage of splitting method
	 * @param s direction aka stage
	 * @param timeStep time step
	 * @param mesh mesh to perform calculation
	 */
	void stage(const int s, const real& timeStep, Mesh& mesh) const {
		for (auto it : mesh) {
			mesh._pdeNew(it) = localGcmStep(
					mesh.matrices(it)->m[s].U1,
					mesh.matrices(it)->m[s].U,
					interpolateValuesAround(mesh, s, it,
							crossingPoints(it, s, timeStep, mesh)));
		}
	}

public: // TODO - private (rewrite test)
	/** Points where characteristics from next time layer cross current time layer */
	PdeVector crossingPoints(const Iterator& it, const int s,
	                         const real& timeStep, const Mesh& mesh) const {
		return -timeStep * linal::diag(mesh.matrices(it)->m[s].L);
	}
	
	
	/**
	 * Interpolate nodal values in specified points.
	 * Interpolated value for k-th point in vector dx
	 * stored in k-th column of returned Matrix.
	 * @param mesh mesh to perform interpolation on
	 * @param s direction
	 * @param it index-iterator of node
	 * @param dx Vector of distances from reference node on which
	 * values should be interpolated
	 * @return Matrix with interpolated nodal values in columns
	 */
	Matrix interpolateValuesAround(const Mesh& mesh, const int s,
	                               const Iterator& it, const PdeVector& dx) const {
		Matrix ans;
		std::vector<PdeVector> src( (size_t) (mesh.borderSize + 1) );
		PdeVector res;
		Iterator shift({0, 0, 0}, mesh.sizes);
		for (int k = 0; k < PdeVector::M; k++) {
			shift(s) = (dx(k) > 0) ? 1 : -1;
			for (int i = 0; i < src.size(); i++) {
				src[(size_t)i] = mesh.pde(it + shift * i);
			}
			EqualDistanceLineInterpolator<PdeVector>::minMaxInterpolate(
					res, src, fabs(dx(k)) / mesh.h(s));
			ans.setColumn(k, res);
		}
		return ans;
	}
	
};


/**
 * Grid-characteristic method specialization for CGAL 2D grid
 */
template<typename TModel, typename TMaterial>
class GridCharacteristicMethod<TModel, Cgal2DGrid, TMaterial> {
public:
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
	 * Do grid-characteristic stage of splitting method
	 * @param s direction aka stage
	 * @param timeStep time step
	 * @param mesh mesh to perform calculation
	 */
	void stage(const int s, const real& timeStep, Mesh& mesh) {
		
		/// calculate spatial derivatives of all mesh pde values ones before stage
		/// in order to use them multiple times while stage calculation 
		estimateGradient(mesh);
		
		/// calculate border nodes
		for (auto borderIter =  mesh.borderBegin(); 
		          borderIter != mesh.borderEnd(); ++borderIter) {
			mesh._pdeNew(*borderIter) = localGcmStep(
					mesh.matrices(*borderIter)->m[s].U1,
					mesh.matrices(*borderIter)->m[s].U,
					interpolateValuesAround(mesh, s, *borderIter,
							crossingPoints(*borderIter, s, timeStep, mesh), true));
			
			borderCorrector(mesh, s, *borderIter);
		}
		
		/// calculate inner nodes
		for (auto innerIter =  mesh.innerBegin(); 
		          innerIter != mesh.innerEnd(); ++innerIter) {
			mesh._pdeNew(*innerIter) = localGcmStep(
					mesh.matrices(*innerIter)->m[s].U1,
					mesh.matrices(*innerIter)->m[s].U,
					interpolateValuesAround(mesh, s, *innerIter,
							crossingPoints(*innerIter, s, timeStep, mesh), false));	
		}
		
	}


private:
	/** Points where characteristics from next time layer cross current time layer */
	PdeVector crossingPoints(const Iterator& it, const int s,
	                         const real& timeStep, const Mesh& mesh) const {
		return -timeStep * linal::diag(mesh.matrices(it)->m[s].L);
	}
	

	/**
	 * Interpolate nodal values in specified points.
	 * Interpolated value for k-th point in vector dx are
	 * stored in k-th column of returned Matrix.
	 * If specified point appears to be out of body
	 * AND it's really border case, matrix column is set to zeros
	 * and outerInvariants is added with the index.
	 * @param mesh mesh to perform interpolation on
	 * @param s direction
	 * @param it index-iterator of node
	 * @param dx Vector of distances from reference node on which
	 * values should be interpolated
	 * @param isBorder is given node border or not
	 * @return Matrix with interpolated nodal values in columns
	 */
	Matrix interpolateValuesAround(const Mesh& mesh, const int s,
	                               const Iterator& it, const PdeVector& dx,
	                               const bool isBorder) {
	    outerInvariants.clear();
		Matrix ans;
		
		for (int k = 0; k < PdeVector::M; k++) {
			
			if (dx(k) == 0) {
			// special for exact hit
				ans.setColumn(k, mesh.pde(it));
				continue;
			}
			
			// point to interpolate respectively to point by given iterator
			Real2 shift({(s == 0) * dx(k), (s == 1) * dx(k)});
			auto t = mesh.findOwnerTriangle(it, shift);
			auto u = PdeVector::Zeros();
			
			if (t.valid) {
			// characteristic hits into body
			// second order interpolate inner value in triangle on current time layer
				u = TriangleInterpolator<PdeVector>::interpolate(
					mesh.coords2d(t.p[0]), mesh.pde(t.p[0]), gradients[mesh.getIndex(t.p[0])],
					mesh.coords2d(t.p[1]), mesh.pde(t.p[1]), gradients[mesh.getIndex(t.p[1])],
					mesh.coords2d(t.p[2]), mesh.pde(t.p[2]), gradients[mesh.getIndex(t.p[2])],
					mesh.coords2d(it) + shift);
					
			} else {
			// characteristic hits out of body
				if (isBorder && linal::dotProduct(mesh.normal(it), shift) > 0) {
				// this is a really border case
				// add outer invariant for border corrector
					outerInvariants.push_back(k);
					u = PdeVector::Zeros();
					
				} else {
				// this means that characteristic hits out of body, 
				// however it is not the border case
					if (!isBorder) {
					// it is from inner node 
						u = interpolateInSpaceTime(mesh, it, shift);
						
					} else {
					// it is from border node that is "inner" on that stage
						u = whenInnerBorderIsOuter(mesh, it, shift);
						
					}
				}
			}
			
			ans.setColumn(k, u);
		}
		
		return ans;
	}
	
	
	/**
	 * Apply border conditions according to Chelnokov's PhD thesis, page 42.
	 * Here used outerInvariants, written at interpolateValuesAround before.
	 */
	void borderCorrector(Mesh& mesh, const int s, const Iterator& it) {
		if (outerInvariants.size() == 0) { return; }
		
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
				i, mesh.matrices(it)->m[s].U1.getColumn(outerInvariants[(size_t)i]));
		}
		
		const auto alpha = linal::solveLinearSystem(B * outerU1, b - B * mesh.pdeNew(it));
		mesh._pdeNew(it) = mesh.pdeNew(it) + outerU1 * alpha;
	}
	
	
	/** 
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border in some point. It's possible either for border and inner nodes.
	 * @note border nodes have to be already calculated
	 */
	PdeVector interpolateInSpaceTime(const Mesh& mesh, 
			const Iterator& it, const Real2& shift) const {
		/// first order interpolate in triangle formed by border points from
		/// current and next time layers (triangle in space-time)

		auto borderEdge = mesh.findCrossingBorder(it, shift);
		Real2 r1 = mesh.coords2d(borderEdge.first);
		Real2 r2 = mesh.coords2d(borderEdge.second);
		Real2 r0 = mesh.coords2d(it);
		Real2 rc = linal::linesIntersection(r1, r2, r0, r0 + shift);
			//< coordinate of border-characteristic intersection
		return TriangleInterpolator<PdeVector>::interpolateInOwner(
			{0, 0}, mesh.pde(borderEdge.first),
			{1, 0}, mesh.pde(borderEdge.second),
			{0, 1}, mesh.pdeNew(borderEdge.first),
			{1, 1}, mesh.pdeNew(borderEdge.second),
			{    linal::length(rc - r1) / linal::length(r2 - r1),
			 1 - linal::length(rc - r0) / linal::length(shift)});
	}
	
	
	/** 
	 * Handle the case when characteristic goes along the line with body border
	 * until body border bends from that line. Then interpolate in space-time 
	 * in that point.
	 * @note border nodes have to be already calculated
	 */
	PdeVector interpolateInSpaceTimeAlongBorder(const Mesh& mesh, 
			const Iterator& it, const Iterator& neighbor, const Real2& shift) const {
		Iterator bend = mesh.findBorderFlexion(it, neighbor);
		Real2 a = mesh.coords2d(it);
		Real2 b = mesh.coords2d(bend);
		real w = linal::length(a - b) / linal::length(shift);
		return (1 - w) * mesh.pdeNew(bend) + w * mesh.pde(bend);
	}
	
	
	/** 
	 * For the case when border node has to be calculated as inner, but
	 * the characteristic hits out of the body
	 */
	PdeVector whenInnerBorderIsOuter(const Mesh& mesh, 
			const Iterator& it, const Real2& shift) const {
		
		auto borderNeighbors = mesh.findBorderNeighbors(it);
		Real2 a = mesh.coords2d(borderNeighbors.first);
		Real2 b = mesh.coords2d(it);
		Real2 c = mesh.coords2d(borderNeighbors.second);
		
		switch (linal::positionRelativeToAngle(a, b, c, b + shift)) {
			case (linal::POSITION::INSIDE):
				return interpolateInSpaceTime(mesh, it, shift);
				break;
			case (linal::POSITION::FIRST_BORDER):
				return interpolateInSpaceTimeAlongBorder(mesh, it, borderNeighbors.first, shift);
				break;
			case (linal::POSITION::SECOND_BORDER):
				return interpolateInSpaceTimeAlongBorder(mesh, it, borderNeighbors.second, shift);
				break;
			case (linal::POSITION::OUTSIDE):
			// really border case is already handled before
				LOG_DEBUG("Bad case: outside the body in whenInnerBorderIsOuter at "
						<< mesh.coords2d(it));
				return mesh.pde(it);
				break;
			default:
				THROW_BAD_METHOD("Unknown position");
		}
	}
	
	
	/** Calculate gradients of mesh pde values at each mesh vertex */
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
			
			gradients[mesh.getIndex(v)] = linal::linearLeastSquares(A, b, W);
		}
	}
	
	
	/** Calculate Hessians of mesh pde values at each mesh vertex */
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
				b(i) = gradients[mesh.getIndex(neighbor)] - gradients[mesh.getIndex(v)];
				
				i++; if (i > 7) { break; }
			}
			
			auto H = linal::linearLeastSquares(A, b); // gradient of gradient
			
			hessians[mesh.getIndex(v)] = {H(0)(0), (H(0)(1) + H(1)(0)) / 2.0, H(1)(1)};
		}
	}
	
		
	/// List of outer Riemann invariants used in borderCorrector.
	/// Invariants are specified by their indices in matrix L.
	std::vector<int> outerInvariants;
	
	/// The storage of gradients of mesh pde values.
	std::vector<PdeGradient> gradients;
	/// The storage of hessians of mesh pde values. (Unused now)
	std::vector<PdeHessian> hessians;
	
	USE_AND_INIT_LOGGER("gcm.Cgal2DGridCharacteristicMethod")
	
};


}

#endif // LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
