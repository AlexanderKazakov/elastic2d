#ifndef LIBGCM_GRIDCHARACTERISTICMETHODCGALGRID_HPP
#define LIBGCM_GRIDCHARACTERISTICMETHODCGALGRID_HPP

#include <lib/mesh/grid/Cgal2DGrid.hpp>
#include <lib/mesh/grid/Cgal3DGrid.hpp> 
// FIXME forward decl?
#include <lib/numeric/gcm/GridCharacteristicMethod.hpp>


namespace gcm {

template<int Dimensionality> struct CgalGrid;

template<> struct CgalGrid<2> { 
	typedef Cgal2DGrid type;
	/// number of neighbor vertices in 2D CGAL grid is unlikely to be more
	/// than 8 (no any guaranty, just estimation)
	static const int MAX_NUMBER_OF_VERTEX_NEIGHBORS = 8;
};

template<> struct CgalGrid<3> {
	typedef Cgal3DGrid type;
	/// number of neighbor vertices in 3D CGAL grid is unlikely to be more
	/// than FIXME (no any guaranty, just estimation)
	static const int MAX_NUMBER_OF_VERTEX_NEIGHBORS = 20;
};


/**
 * Grid-characteristic method implementation for
 * meshes based on CGAL 2D and 3D grids
 */
template<typename TModel, typename TMaterial, int Dimensionality>
class GridCharacteristicMethodCgalGrid {
public:
	typedef typename CgalGrid<Dimensionality>::type            Grid;
	typedef TModel                                             Model;
	typedef DefaultMesh<Model, Grid, TMaterial>                Mesh;
	typedef typename Mesh::Matrix                              Matrix;
	typedef typename Mesh::PdeVector                           PdeVector;
	typedef typename Mesh::Iterator                            Iterator;
	typedef typename Mesh::BORDER_CONDITION                    BORDER_CONDITION;
	typedef typename Mesh::Cell                                Cell;
	
	typedef linal::VECTOR<Dimensionality, PdeVector>           PdeGradient;
	typedef linal::SYMMETRIC_MATRIX<Dimensionality, PdeVector> PdeHessian;
	typedef linal::Vector<Dimensionality>                      RealD;
	
	static const int OUTER_NUMBER = BORDER_CONDITION::OUTER_NUMBER;
	static const int PDE_SIZE = Model::PDE_SIZE;
	typedef linal::Matrix<PDE_SIZE, OUTER_NUMBER> OuterU1Matrix;
	
	static const int MAX_NUMBER_OF_VERTEX_NEIGHBORS = 
			CgalGrid<Dimensionality>::MAX_NUMBER_OF_VERTEX_NEIGHBORS;

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
			RealD shift = RealD::Zeros();
			shift(s) = dx(k);
			Cell t = mesh.findOwnerCell(it, shift);
			PdeVector u = PdeVector::Zeros();
			
			if (t.valid) {
			// characteristic hits into body
			// second order interpolate inner value in triangle on current time layer
				u = interpolateInSpace(mesh, mesh.coordsD(it) + shift, t);
					
			} else {
			// characteristic hits out of body
				if (isBorder && linal::dotProduct(mesh.normal(it), shift) > 0) {
				// this is a really border case
				// add outer invariant for border corrector
//					outerInvariants.push_back(k);
					u = PdeVector::Zeros();
					
				} else {
				// this means that characteristic hits out of body, 
				// however it is not the border case
					if (!isBorder) {
					// it is from inner node 
//						u = interpolateInSpaceTime(mesh, it, shift);
						
					} else {
					// it is from border node that is "inner" on that stage
//						u = whenInnerBorderIsOuter(mesh, it, shift);
						
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
				<< "in border node at " << mesh.coordsD(it));
			return;
		}
		
		const BORDER_CONDITION* borderCondition = mesh.getBorderCondition(it);
		if (borderCondition == nullptr) { return; } // non-reflection
		
		const RealD normal = mesh.normal(it);
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
	
	
	/** Interpolate PdeVector from space on current time layer (2D case) */
	PdeVector interpolateInSpace(const Mesh& mesh, const Real2& query, const Cell& t) const {
		return TriangleInterpolator<PdeVector>::interpolate(
				mesh.coordsD(t(0)), mesh.pde(t(0)), gradients[mesh.getIndex(t(0))],
				mesh.coordsD(t(1)), mesh.pde(t(1)), gradients[mesh.getIndex(t(1))],
				mesh.coordsD(t(2)), mesh.pde(t(2)), gradients[mesh.getIndex(t(2))],
				query);
	}
	
	
	/** Interpolate PdeVector from space on current time layer (3D case) */
	PdeVector interpolateInSpace(const Mesh& mesh, const Real3& query, const Cell& t) const {
		return TetrahedronInterpolator<PdeVector>::interpolate(
				mesh.coordsD(t(0)), mesh.pde(t(0)), gradients[mesh.getIndex(t(0))],
				mesh.coordsD(t(1)), mesh.pde(t(1)), gradients[mesh.getIndex(t(1))],
				mesh.coordsD(t(2)), mesh.pde(t(2)), gradients[mesh.getIndex(t(2))],
				mesh.coordsD(t(3)), mesh.pde(t(3)), gradients[mesh.getIndex(t(3))],
				query);
	}
	
	
	/** 
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border in some point. It's possible either for border and inner nodes.
	 * @note border nodes must be already calculated
	 */
	PdeVector interpolateInSpaceTime(const Mesh& mesh, 
			const Iterator& it, const Real2& shift) const {
		/// 2D case
		/// first order interpolate in triangle formed by border points from
		/// current and next time layers (triangle in space-time)

		auto borderEdge = mesh.findCrossingBorder(it, shift);
		Real2 r1 = mesh.coordsD(borderEdge.first);
		Real2 r2 = mesh.coordsD(borderEdge.second);
		Real2 r0 = mesh.coordsD(it);
		Real2 rc = linal::linesIntersection(r1, r2, r0, r0 + shift);
				//< coordinate of border-characteristic intersection
		
		return TriangleInterpolator<PdeVector>::interpolateInOwner(
				// current time layer
				{0, 0}, mesh.pde(borderEdge.first),
				{1, 0}, mesh.pde(borderEdge.second),
				// next time layer
				{0, 1}, mesh.pdeNew(borderEdge.first),
				{1, 1}, mesh.pdeNew(borderEdge.second),
				// query in space-time
				{    linal::length(rc - r1) / linal::length(r2 - r1),
				 1 - linal::length(rc - r0) / linal::length(shift)});
	}
	
	
	/** 
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border in some point. It's possible either for border and inner nodes.
	 * @note border nodes must be already calculated
	 */
	PdeVector interpolateInSpaceTime(const Mesh& mesh, 
			const Iterator& it, const Real3& shift) const {
		/// 3D case
		/// first order interpolate in tetrahedron formed by border points from
		/// current and next time layers (tetrahedron in space-time)

		auto borderFace = mesh.findCrossingBorder(it, shift);
		Real3 r1 = mesh.coordsD(borderFace(0));
		Real3 r2 = mesh.coordsD(borderFace(1));
		Real3 r3 = mesh.coordsD(borderFace(2));
		Real3 r0 = mesh.coordsD(it);
		Real3 rc = linal::lineWithFlatIntersection(r1, r2, r3, r0, r0 + shift);
				//< coordinate of border-characteristic intersection
		
		return TetrahedronInterpolator<PdeVector>::interpolateInOwner(
				// current time layer
				{0, 0, 0}, mesh.pde(borderFace(0)),
				{1, 0, 0}, mesh.pde(borderFace(1)),
				{0, 1, 0}, mesh.pde(borderFace(2)),
				// next time layer
				{0, 0, 1}, mesh.pdeNew(borderFace(0)),
				{1, 0, 1}, mesh.pdeNew(borderFace(1)),
				{0, 1, 1}, mesh.pdeNew(borderFace(2)),
				// query in space-time
				{    linal::length(rc - r1) / linal::length(r2 - r1),
					 linal::length(rc - r1) / linal::length(r3 - r1),
				 1 - linal::length(rc - r0) / linal::length(shift)});
	}
	
	
	/** 
	 * For the case when border node has to be calculated as inner, but
	 * the characteristic hits out of the body (2D case)
	 */
	PdeVector whenInnerBorderIsOuter(const Mesh& mesh, 
			const Iterator& it, const Real2& shift) const {
		
		auto borderNeighbors = mesh.findBorderNeighbors(it);
		Real2 a = mesh.coordsD(borderNeighbors.first);
		Real2 b = mesh.coordsD(it);
		Real2 c = mesh.coordsD(borderNeighbors.second);
		
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
						<< mesh.coordsD(it));
				return mesh.pde(it);
				break;
			default:
				THROW_BAD_METHOD("Unknown position");
		}
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
		Real2 a = mesh.coordsD(it);
		Real2 b = mesh.coordsD(bend);
		real w = linal::length(a - b) / linal::length(shift);
		return (1 - w) * mesh.pdeNew(bend) + w * mesh.pde(bend);
	}
	
	
	/** 
	 * For the case when border node has to be calculated as inner, but
	 * the characteristic hits out of the body (3D case)
	 */
	PdeVector whenInnerBorderIsOuter(const Mesh& mesh, 
			const Iterator& it, const Real3& /*shift*/) const {	
		LOG_DEBUG("TODO - whenInnerBorderIsOuter at " << mesh.coordsD(it));
		return mesh.pde(it);
	}
	
	
	// FIXME - move to separate class
	/** Calculate gradients of mesh pde values at each mesh vertex */
	void estimateGradient(const Mesh& mesh) {
		gradients.resize(mesh.sizeOfAllNodes());

		// SLE matrix
		linal::Matrix<MAX_NUMBER_OF_VERTEX_NEIGHBORS, Dimensionality> A;
		// SLE right part
		linal::VECTOR<MAX_NUMBER_OF_VERTEX_NEIGHBORS, PdeVector> b;
		// weight matrix for linear least squares method
		linal::DiagonalMatrix<MAX_NUMBER_OF_VERTEX_NEIGHBORS> W;
		
		for (const auto v : mesh) {
			linal::clear(A); linal::clear(b); linal::clear(W);
			
			/// estimate gradient from its definition:
			/// \f$  (\nabla f, \vec{a} - \vec{b}) = f(a) - f(b),  \f$ 
			/// applying it to each neighbor of the vertex
			const auto neighbors = mesh.findNeighborVertices(v);
			int i = 0;
			for (const auto neighbor : neighbors) {
				RealD d = mesh.coordsD(neighbor) - mesh.coordsD(v);
				A.setRow(i, d);
				W(i) = linal::length(d);
				b(i) = mesh.pde(neighbor) - mesh.pde(v);
				
				i++; if (i == MAX_NUMBER_OF_VERTEX_NEIGHBORS) { break; }
			}
			
			gradients[mesh.getIndex(v)] = linal::linearLeastSquares(A, b, W);
		}
	}
	
	
	/** Calculate Hessians of mesh pde values at each mesh vertex */
	void estimateHessian(const Mesh& mesh) {
	/// @note gradients must be already calculated
	/// the function is currently unused because third order is not 
	/// so straightforward as second order
		hessians.resize(mesh.sizeOfAllNodes());
		assert_eq(hessians.size(), gradients.size());

		// SLE matrix
		linal::Matrix<MAX_NUMBER_OF_VERTEX_NEIGHBORS, Dimensionality> A;
		// SLE right part
		linal::VECTOR<MAX_NUMBER_OF_VERTEX_NEIGHBORS, PdeGradient> b;
		// weight matrix for linear least squares method
		linal::DiagonalMatrix<MAX_NUMBER_OF_VERTEX_NEIGHBORS> W;

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
				W(i) = linal::length(d);
				b(i) = gradients[mesh.getIndex(neighbor)] - gradients[mesh.getIndex(v)];
				
				i++; if (i == MAX_NUMBER_OF_VERTEX_NEIGHBORS) { break; }
			}
			
			auto H = linal::linearLeastSquares(A, b, W); // gradient of gradient
			
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
	
	USE_AND_INIT_LOGGER("gcm.GridCharacteristicMethodCgalGrid")
	
};


/**
 * Grid-characteristic method for meshes based on Cgal2DGrid
 */
template<typename TModel, typename TMaterial>
class GridCharacteristicMethod<TModel, Cgal2DGrid, TMaterial> :
		public GridCharacteristicMethodCgalGrid<TModel, TMaterial, 2> { };
/**
 * Grid-characteristic method for meshes based on Cgal3DGrid
 */
template<typename TModel, typename TMaterial>
class GridCharacteristicMethod<TModel, Cgal3DGrid, TMaterial> :
		public GridCharacteristicMethodCgalGrid<TModel, TMaterial, 3> { };


}

#endif // LIBGCM_GRIDCHARACTERISTICMETHODCGALGRID_HPP
