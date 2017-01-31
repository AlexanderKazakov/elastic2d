#ifndef LIBGCM_SIMPLEX_GRIDCHARACTERISTICMETHOD_HPP
#define LIBGCM_SIMPLEX_GRIDCHARACTERISTICMETHOD_HPP

#include <libgcm/util/infrastructure/infrastructure.hpp>
#include <libgcm/grid/AbstractGrid.hpp>
#include <libgcm/util/math/Differentiation.hpp>
#include <libgcm/util/math/interpolation/interpolation.hpp>
#include <libgcm/util/math/GridCharacteristicMethod.hpp>


namespace gcm {
namespace simplex {


class GridCharacteristicMethodBase {
public:
	virtual void beforeStage(
			const AbstractGrid& mesh_) = 0;
	virtual void borderStage(
			const int s, const real timeStep, AbstractGrid& mesh_) = 0;
	virtual void innerStage(
			const int s, const real timeStep, AbstractGrid& mesh_) = 0;
};


/**
 * Grid-characteristic method for meshes based on SimplexGrid
 */
template<typename Mesh>
class GridCharacteristicMethod : public GridCharacteristicMethodBase {
public:
	typedef typename Mesh::Matrix                              Matrix;
	typedef typename Mesh::PdeVector                           PdeVector;
	typedef typename Mesh::Iterator                            Iterator;
	typedef typename Mesh::Cell                                Cell;
	
	typedef Differentiation<Mesh>                              DIFFERENTIATION;
	typedef typename DIFFERENTIATION::PdeGradient              PdeGradient;
	typedef typename DIFFERENTIATION::PdeHessian               PdeHessian;
	typedef typename DIFFERENTIATION::RealD                    RealD;
	
	typedef typename Mesh::Model                               Model;
	static const int OUTER_NUMBER = Model::OUTER_NUMBER;
	
	
	virtual void beforeStage(
			const AbstractGrid& mesh_) override {
		/// calculate spatial derivatives of all mesh pde values ones before stage
		/// in order to use them multiple times while stage calculation 
		const Mesh& mesh = dynamic_cast<const Mesh&>(mesh_);
		DIFFERENTIATION::estimateGradient(mesh, gradients);
	}
	
	
	/**
	 * Do grid-characteristic stage of splitting method on contact and border nodes
	 */
	virtual void borderStage(
			const int s, const real timeStep, AbstractGrid& mesh_) override {
		Mesh& mesh = dynamic_cast<Mesh&>(mesh_);
		RealD direction = mesh.getCalculationBasis().getColumn(s);
		
		/// calculate inner waves of contact nodes in Riemann variables
		/// then, it will be handled by ContactCorrector
		for (auto contactIter = mesh.contactBegin();
		          contactIter < mesh.contactEnd(); ++contactIter) {
			mesh._pdeNew(*contactIter) = linal::diagonalMultiply(
					mesh.matrices(*contactIter)->m[s].U,
					interpolateValuesAroundBorderNode(
							mesh, direction, *contactIter,
							crossingPoints(*contactIter, s, timeStep, mesh)));
			mesh._waveIndices(*contactIter) = innerWaves;
		}
		
		
		/// calculate inner waves of border nodes in Riemann variables
		/// then, it will be handled by BorderCorrector
		for (auto borderIter = mesh.borderBegin();
		          borderIter < mesh.borderEnd(); ++borderIter) {
			mesh._pdeNew(*borderIter) = linal::diagonalMultiply(
					mesh.matrices(*borderIter)->m[s].U,
					interpolateValuesAroundBorderNode(
							mesh, direction, *borderIter,
							crossingPoints(*borderIter, s, timeStep, mesh)));
			mesh._waveIndices(*borderIter) = innerWaves;
		}
	}
	
	
	/**
	 * Do grid-characteristic stage of splitting method on inner nodes
	 * @note contact and border nodes must be already calculated
	 * @param s number of stage (GcmMatrix number)
	 * @param timeStep time step
	 * @param mesh mesh to perform calculation
	 */
	virtual void innerStage(
			const int s, const real timeStep, AbstractGrid& mesh_) override {
		Mesh& mesh = dynamic_cast<Mesh&>(mesh_);
		RealD direction = mesh.getCalculationBasis().getColumn(s);
		
		/// calculate inner nodes in PDE variables
		for (auto innerIter = mesh.innerBegin(); 
		          innerIter < mesh.innerEnd(); ++innerIter) {
			mesh._pdeNew(*innerIter) = localGcmStep(
					mesh.matrices(*innerIter)->m[s].U1,
					mesh.matrices(*innerIter)->m[s].U,
					interpolateValuesAroundInnerNode(
							mesh, direction, *innerIter,
							crossingPoints(*innerIter, s, timeStep, mesh)));
		}
	}
	
	
private:
	/** Points where characteristics from next time layer cross current time layer */
	PdeVector crossingPoints(const Iterator& it, const int s,
	                         const real timeStep, const Mesh& mesh) const {
		return -timeStep * linal::diag(mesh.matrices(it)->m[s].L);
	}
	
	
	/**
	 * Interpolate PDE vectors in specified points.
	 * Interpolated vector for k-th point in vector dx are
	 * stored in k-th column in returned Matrix.
	 * If specified point appears to be out of body,
	 * which is possible for border nodes, matrix column is set to zeros.
	 * @param mesh mesh to perform interpolation on
	 * @param direction direction of line to find values along
	 * @param it index-iterator of node
	 * @param dx Vector of distances from the node on which
	 * values should be interpolated
	 * @return Matrix with interpolated nodal values in columns
	 */
	Matrix interpolateValuesAroundBorderNode(
			const Mesh& mesh, const RealD direction,
			const Iterator& it, const PdeVector& dx) {
		innerWaves.clear();
		Matrix ans = Matrix::Zeros();
		
		for (int k = 0; k < PdeVector::M; k++) {
			
			if (dx(k) == 0) {
			// special for exact hit
				ans.setColumn(k, mesh.pde(it));
				continue;
			}
			
			// point to interpolate respectively to point by given iterator
			RealD shift = direction * dx(k);
			Cell t = mesh.findCellCrossedByTheRay(it, shift);
			if (t.n == t.N) {
			// characteristic hits into the body
				ans.setColumn(k,
						interpolateInSpace(mesh, mesh.coordsD(it) + shift, t));
				innerWaves.push_back(k);
			}
			
		}
		return ans;
	}
	
	
	/**
	 * @see interpolateValuesAroundBorderNode
	 * Calculation of inner nodes differs from border nodes in two points:
	 * 1. Interpolation can be performed not in space only, but also 
	 *    in space-time using values from border nodes on the next time layer.
	 * 2. Thus, there can not be outer characteristic hits from inner nodes.
	 */
	Matrix interpolateValuesAroundInnerNode(
			const Mesh& mesh, const RealD direction,
			const Iterator& it, const PdeVector& dx) {
		Matrix ans = Matrix::Zeros();
		
		for (int k = 0; k < PdeVector::M; k++)  {
			
			if (dx(k) == 0) {
			// special for exact hit
				ans.setColumn(k, mesh.pde(it));
				continue;
			}
			
			// point to interpolate respectively to point by given iterator
			RealD shift = direction * dx(k);
			Cell t = mesh.findCellCrossedByTheRay(it, shift);
			
			PdeVector u = PdeVector::Zeros();
			if (t.n == t.N) {
			// characteristic hits into the body
				u = interpolateInSpace(mesh, mesh.coordsD(it) + shift, t);
			} else if (t.n == t.N - 1) {
			// characteristic hits out of body going throughout border face
				u = interpolateInSpaceTime(mesh, it, shift, t);
			} else if (t.n == t.N - 2) {
			// exact hit to border edge (3D) or point (2D)
				u = interpolateInSpaceTime1D(mesh, it, shift, t);
			} else {
				THROW_BAD_METHOD("Outer hit from inner node?");
			}
			ans.setColumn(k, u);
			
		}
		return ans;
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
		/// coordinate of border-characteristic intersection
		Real3 rc = linal::lineWithFlatIntersection(r1, r2, r3, r0, r0 + shift);
		/// find local coordinates of (rc - r1) in basis {e1, e2}
		Real3 e1 = r2 - r1, e2 = r3 - r1, a = rc - r1;
		linal::Matrix<3, 2> A = {
				e1(0), e2(0),
				e1(1), e2(1),
				e1(2), e2(2),
		};
		Real2 w = linal::linearLeastSquares(A, a);
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
				{w(0), w(1), 1 - linal::length(rc - r0) / linal::length(shift)});
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
		THROW_UNSUPPORTED("TODO");
	}
	
	
	/// List of inner Riemann invariants after border/contact node calculation.
	/// Used by border/contact correctors afterwards
	typename Mesh::WaveIndices innerWaves;
	
	/// The storage of gradients of mesh pde values.
	std::vector<PdeGradient> gradients;
	
	USE_AND_INIT_LOGGER("gcm.simplex.GridCharacteristicMethod")
};


} // namespace simplex
} // namespace gcm


#endif // LIBGCM_SIMPLEX_GRIDCHARACTERISTICMETHOD_HPP
