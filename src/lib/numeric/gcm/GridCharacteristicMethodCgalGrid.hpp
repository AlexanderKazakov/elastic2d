#ifndef LIBGCM_GRIDCHARACTERISTICMETHODCGALGRID_HPP
#define LIBGCM_GRIDCHARACTERISTICMETHODCGALGRID_HPP

#include <lib/mesh/grid/Cgal2DGrid.hpp>
#include <lib/mesh/grid/Cgal3DGrid.hpp> 
#include <lib/numeric/gcm/GridCharacteristicMethod.hpp>
#include <lib/numeric/gcm/Differentiation.hpp>


namespace gcm {

template<int Dimensionality> struct CgalGrid;

template<> struct CgalGrid<2> { 
	typedef Cgal2DGrid type;
};

template<> struct CgalGrid<3> {
	typedef Cgal3DGrid type;
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
	
	typedef Differentiation<Mesh>                              DIFFERENTIATION;
	typedef linal::VECTOR<Dimensionality, PdeVector>           PdeGradient;
	typedef linal::SYMMETRIC_MATRIX<Dimensionality, PdeVector> PdeHessian;
	typedef linal::Vector<Dimensionality>                      RealD;
	
	static const int OUTER_NUMBER = BORDER_CONDITION::OUTER_NUMBER;
	static const int PDE_SIZE = Model::PDE_SIZE;
	typedef linal::Matrix<PDE_SIZE, OUTER_NUMBER> OuterU1Matrix;
	
	/**
	 * Do grid-characteristic stage of splitting method
	 * @param s number of stage (GcmMatrix number)
	 * @param timeStep time step
	 * @param mesh mesh to perform calculation
	 * @param direction direction of line to perform stage along
	 */
	void stage(const int s, const real timeStep, Mesh& mesh, const RealD direction) {
		assert_eq(linal::length(direction), 1);
		
		/// calculate spatial derivatives of all mesh pde values ones before stage
		/// in order to use them multiple times while stage calculation 
		DIFFERENTIATION::estimateGradient(mesh, gradients);
		
		/// calculate border nodes
		LOG_INFO("Start calculate border nodes");
//		size_t counter = 0;
//		#pragma omp parallel for
		for (auto borderIter = mesh.borderBegin(); 
		          borderIter < mesh.borderEnd(); ++borderIter) {
			mesh._pdeNew(*borderIter) = localGcmStep(
					mesh.matrices(*borderIter)->m[s].U1,
					mesh.matrices(*borderIter)->m[s].U,
					interpolateValuesAround(mesh, direction, *borderIter,
							crossingPoints(*borderIter, s, timeStep, mesh), true));
			
			borderCorrector(mesh, s, direction, *borderIter);
			
//			if (++counter % 20000 == 0) {
//				LOG_INFO(counter << " border nodes have been calculated");
//			}
		}
		
		/// calculate inner nodes
		LOG_INFO("Start calculate inner nodes");
//		counter = 0;
//		#pragma omp parallel for
		for (auto innerIter = mesh.innerBegin(); 
		          innerIter < mesh.innerEnd(); ++innerIter) {
			mesh._pdeNew(*innerIter) = localGcmStep(
					mesh.matrices(*innerIter)->m[s].U1,
					mesh.matrices(*innerIter)->m[s].U,
					interpolateValuesAround(mesh, direction, *innerIter,
							crossingPoints(*innerIter, s, timeStep, mesh), false));
			
//			if (++counter % 20000 == 0) {
//				LOG_INFO(counter << " inner nodes have been calculated");
//			}
		}
		
	}


private:
	/** Points where characteristics from next time layer cross current time layer */
	PdeVector crossingPoints(const Iterator& it, const int s,
	                         const real timeStep, const Mesh& mesh) const {
		return -timeStep * linal::diag(mesh.matrices(it)->m[s].L);
	}
	

	/**
	 * Interpolate nodal values in specified points.
	 * Interpolated value for k-th point in vector dx are
	 * stored in k-th column of returned Matrix.
	 * If specified point appears to be out of body
	 * AND it is really border case, matrix column is set to zeros
	 * and outerInvariants is added with the index.
	 * @param mesh mesh to perform interpolation on
	 * @param direction direction of line to find values along
	 * @param it index-iterator of node
	 * @param dx Vector of distances from reference node on which
	 * values should be interpolated
	 * @param isBorder is given node border or not
	 * @return Matrix with interpolated nodal values in columns
	 */
	Matrix interpolateValuesAround(const Mesh& mesh, const RealD direction,
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
			RealD shift = direction * dx(k);
			Cell t = mesh.findOwnerCell(it, shift);
			PdeVector u = PdeVector::Zeros();
			
			if (t.n == t.N) {
			// characteristic hits into body
			// second order interpolate inner value in triangle on current time layer
				u = interpolateInSpace(mesh, mesh.coordsD(it) + shift, t);
					
			} else if (t.n == t.N - 1) {
			// characteristic hits out of body going throughout border face
				if (!isBorder) { // TODO
				// if node is border, the base of interpolation can be not calculated yet
					u = interpolateInSpaceTime(mesh, it, shift, t);
				}
			
			} else if (t.n == t.N - 2) {
			// characteristic hits out of body going throughout corner of border face
				if (!isBorder) { // TODO
				// if node is border, the base of interpolation can be not calculated yet
					u = interpolateInSpaceTime1D(mesh, it, shift, t);
				}
			
			} else {
				/*if (!isBorder) {
					LOG_INFO("Missed hit from inner node at " << mesh.coordsD(it)
							<< "t.n == " << t.n);
				}*/ // FIXME
				/*assert_true(isBorder);*/
			// now, all inner and part of border cases are calculated	
				if (isBorder /*&&
					linal::dotProduct(mesh.normal(it),
							linal::normalize(shift)) > cos(M_PI / 4 + EQUALITY_TOLERANCE)*/) {
				// this is really border case
				// add outer invariant for border corrector
					outerInvariants.push_back(k);
					u = PdeVector::Zeros();
				
				} else {
				// it can not be calculated as border, because it would be 
				// numerically unstable (smth like incompatible border conditions)
//					u = mesh.pde(it);
					/*LOG_DEBUG("Bad case at " << mesh.coordsD(it) << std::endl
							<< "shift: " << shift << std::endl);*/
					
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
	void borderCorrector(Mesh& mesh, const int /*s*/, const RealD direction,
			const Iterator& it) {
		
		if (outerInvariants.size() == 0) { return; }
		
		if (outerInvariants.size() == 2 * OUTER_NUMBER) {
		/// no space - no waves
			mesh._pdeNew(it) = mesh.pde(it);
			return;
		}
		
		if (outerInvariants.size() != OUTER_NUMBER) {
//			LOG_INFO("Bad case: " << outerInvariants.size() << " outer invariants "
//					<< "in border node at " << mesh.coordsD(it));
			return;
		}
		
		const BORDER_CONDITION* borderCondition = mesh.getBorderCondition(it);
		if (borderCondition == nullptr) { return; } // non-reflection
		
		const RealD normal = mesh.normal(it);
		const auto B = borderCondition->B(normal);
		const auto b = borderCondition->b();
		
		Matrix u1AlongBorderNormal;
		Model::constructEigenvectors(u1AlongBorderNormal,
				mesh.material(it), linal::createLocalBasis(
						normal * Utils::sign(linal::dotProduct(normal, direction))));
		
		OuterU1Matrix outerU1;		 
		for (int i = 0; i < OUTER_NUMBER; i++) {
			outerU1.setColumn(i, 
					u1AlongBorderNormal.getColumn(outerInvariants[(size_t)i]));
		}
		
		const auto alpha = linal::solveLinearSystem(B * outerU1, b - B * mesh.pdeNew(it));
		mesh._pdeNew(it) = mesh.pdeNew(it) + outerU1 * alpha;
	}
	
	
	/** Interpolate PdeVector from space on current time layer (2D case) */
	PdeVector interpolateInSpace(const Mesh& mesh, const Real2& query, const Cell& c) const {
		return TriangleInterpolator<PdeVector>::interpolate(
				mesh.coordsD(c(0)), mesh.pde(c(0)), gradients[mesh.getIndex(c(0))],
				mesh.coordsD(c(1)), mesh.pde(c(1)), gradients[mesh.getIndex(c(1))],
				mesh.coordsD(c(2)), mesh.pde(c(2)), gradients[mesh.getIndex(c(2))],
				query);
	}
	
	
	/** Interpolate PdeVector from space on current time layer (3D case) */
	PdeVector interpolateInSpace(const Mesh& mesh, const Real3& query, const Cell& c) const {
		return TetrahedronInterpolator<PdeVector>::interpolate(
				mesh.coordsD(c(0)), mesh.pde(c(0)), gradients[mesh.getIndex(c(0))],
				mesh.coordsD(c(1)), mesh.pde(c(1)), gradients[mesh.getIndex(c(1))],
				mesh.coordsD(c(2)), mesh.pde(c(2)), gradients[mesh.getIndex(c(2))],
				mesh.coordsD(c(3)), mesh.pde(c(3)), gradients[mesh.getIndex(c(3))],
				query);
	}
	
	
	/** 
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border in some point. It's possible either for border and inner nodes.
	 * @note border nodes must be already calculated
	 */
	PdeVector interpolateInSpaceTime(const Mesh& mesh, 
			const Iterator& it, const Real2& shift, const Cell& borderEdge) const {
		/// 2D case
		/// first order interpolate in triangle formed by border points from
		/// current and next time layers (triangle in space-time)

		Real2 r1 = mesh.coordsD(borderEdge(0));
		Real2 r2 = mesh.coordsD(borderEdge(1));
		Real2 r0 = mesh.coordsD(it);
		Real2 rc = linal::linesIntersection(r1, r2, r0, r0 + shift);
				//< coordinate of border-characteristic intersection
		
		return TriangleInterpolator<PdeVector>::interpolateInOwner(
				// current time layer
				{0, 0}, mesh.pde(borderEdge(0)),
				{1, 0}, mesh.pde(borderEdge(1)),
				// next time layer
				{0, 1}, mesh.pdeNew(borderEdge(0)),
				{1, 1}, mesh.pdeNew(borderEdge(1)),
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
			const Iterator& it, const Real3& shift, const Cell& borderFace) const {
		/// 3D case
		/// first order interpolate in tetrahedron formed by border points from
		/// current and next time layers (tetrahedron in space-time)

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
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border exactly through the border vertex. It's possible for 2D case.
	 * It's possible either for border and inner nodes.
	 * @note border nodes must be already calculated
	 */
	PdeVector interpolateInSpaceTime1D(const Mesh& mesh, 
			const Iterator& it, const Real2& shift, const Cell& borderVertex) const {
		/// 2D case
		/// first order interpolate in the line formed by crossed point
		/// at current and next time layers (line in space-time)

		auto bv = borderVertex(0);
		Real2 rv = mesh.coordsD(bv);
		Real2 r0 = mesh.coordsD(it);
		real w = linal::length(rv - r0) / linal::length(shift);
		return mesh.pde(bv) * w + mesh.pdeNew(it) * (1 - w);
	}
	
	
	/** 
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border exactly through the border vertex. It's possible for 2D case.
	 * It's possible either for border and inner nodes.
	 * @note border nodes must be already calculated
	 */
	PdeVector interpolateInSpaceTime1D(const Mesh& /*mesh*/, 
			const Iterator& /*it*/, const Real3& /*shift*/, const Cell& /*borderEdge*/) const {
		return PdeVector::Zeros(); // FIXME
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
