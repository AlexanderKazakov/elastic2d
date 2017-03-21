#ifndef LIBGCM_SIMPLEX_GRIDCHARACTERISTICMETHODINRIEMANNINVARIANTS_HPP
#define LIBGCM_SIMPLEX_GRIDCHARACTERISTICMETHODINRIEMANNINVARIANTS_HPP

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
			const int s, AbstractGrid& mesh_) = 0;
	virtual void contactAndBorderStage(
			const int s, const real timeStep, AbstractGrid& mesh_) = 0;
	virtual void innerStage(
			const int s, const real timeStep, AbstractGrid& mesh_) = 0;
	virtual void afterStage(
			const int s, AbstractGrid& mesh_) = 0;
};


/**
 * Grid-characteristic method for meshes based on SimplexGrid.
 * The approach is to calculate and advect along characteristics
 * scalar Riemann-invariants not PDE-vectors
 * @see GridCharacteristicMethodInPdeVectors -- an opposite approach
 */
template<typename Mesh>
class GridCharacteristicMethodInRiemannInvariants :
		public GridCharacteristicMethodBase {
public:
	typedef typename Mesh::Matrix                              Matrix;
	typedef typename Mesh::PdeVector                           PdeVector;
	typedef typename Mesh::PdeVariables                        PdeVariables;
	typedef typename Mesh::GCM_MATRICES                        GcmMatrices;
	typedef typename Mesh::Iterator                            Iterator;
	typedef typename Mesh::Cell                                Cell;
	typedef typename Mesh::WaveIndices                         WaveIndices;
	
	typedef Differentiation<Mesh>                              DIFFERENTIATION;
	typedef typename DIFFERENTIATION::PdeGradient              PdeGradient;
	typedef typename DIFFERENTIATION::PdeHessian               PdeHessian;
	typedef typename DIFFERENTIATION::RealD                    RealD;
	
	typedef typename Mesh::Model                               Model;
	static const int OUTER_NUMBER = Model::OUTER_NUMBER;
	static const int DIMENSIONALITY = Mesh::DIMENSIONALITY;
	static const int CELL_POINTS_NUMBER = Mesh::CELL_POINTS_NUMBER;
	
	typedef real	                              RiemannInvariant;
	typedef linal::VECTOR<
			DIMENSIONALITY, RiemannInvariant>	 RiemannInvariantGradient;
	
	
	virtual void beforeStage(
			const int s, AbstractGrid& mesh_) override {
		Mesh& mesh = dynamic_cast<Mesh&>(mesh_);
		savedPdeTimeLayer = mesh.getPdeVariablesStorage();
		/// switch to Riemann invariants for that stage
		for (const Iterator it : mesh) {
			mesh._pde(it) = (*mesh.matrices(it))(s).U * mesh.pde(it);
		}
		/// calculate spatial derivatives of all mesh Riemann invariants ones before stage
		/// in order to use them multiple times while stage calculation 
		DIFFERENTIATION::estimateGradient(mesh, gradients);
	}
	
	
	/**
	 * Do grid-characteristic stage of splitting method on contact and border nodes
	 */
	virtual void contactAndBorderStage(
			const int s, const real timeStep, AbstractGrid& mesh_) override {
		Mesh& mesh = dynamic_cast<Mesh&>(mesh_);
		
		/// calculate inner waves of contact nodes
		for (auto iter = mesh.contactBegin(); iter < mesh.contactEnd(); ++iter) {
			const GcmMatrices& gcmMatrices = *mesh.matrices(*iter);
			const RealD direction = gcmMatrices.basis.getColumn(s);
			mesh._pdeNew(s, *iter) = interpolateValuesAround(
					s, mesh, direction, *iter,
					crossingPoints(*iter, s, timeStep, mesh), false);
			mesh._waveIndices(*iter) = outerInvariants;
		}
		
		/// calculate inner waves of border nodes
		for (auto iter = mesh.borderBegin(); iter < mesh.borderEnd(); ++iter) {
			const GcmMatrices& gcmMatrices = *mesh.matrices(*iter);
			const RealD direction = gcmMatrices.basis.getColumn(s);
			mesh._pdeNew(s, *iter) = interpolateValuesAround(
					s, mesh, direction, *iter,
					crossingPoints(*iter, s, timeStep, mesh), false);
			mesh._waveIndices(*iter) = outerInvariants;
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
		const RealD direction = mesh.getInnerCalculationBasis().getColumn(s);
		
		/// calculate inner nodes
		for (auto iter = mesh.innerBegin(); iter < mesh.innerEnd(); ++iter) {
			mesh._pdeNew(s, *iter) = interpolateValuesAround(
					s, mesh, direction, *iter,
					crossingPoints(*iter, s, timeStep, mesh), true);
			assert_eq(outerInvariants.size(), 0);
		}
	}
	
	
	virtual void afterStage(
			const int s, AbstractGrid& mesh_) override {
		Mesh& mesh = dynamic_cast<Mesh&>(mesh_);
		/// set PDE-vectors on old time layer to its values saved before the stage
		mesh.swapPdeVariablesStorage(savedPdeTimeLayer);
		/// switch back from Riemann invariants to PDE-variables
		for (const Iterator it : mesh) {
			mesh._pdeNew(s, it) = (*mesh.matrices(it))(s).U1 * mesh.pdeNew(s, it);
		}
	}
	
	
private:
	/** Points where characteristics from next time layer cross current time layer */
	PdeVector crossingPoints(const Iterator& it, const int s,
	                         const real timeStep, const Mesh& mesh) const {
		return -timeStep * linal::diag(mesh.matrices(it)->m[s].L);
	}
	
	
	/**
	 * Interpolate Riemann invariants in specified points.
	 * If specified point appears to be out of body
	 * AND it is really border case, invariant is set to zero
	 * and outerInvariants is added with the index.
	 * @param s stage
	 * @param mesh mesh to perform interpolation on
	 * @param direction direction of line to find values along
	 * @param it index-iterator of node
	 * @param dx Vector of distances from reference node on which
	 * values should be interpolated
	 * @param canInterpolateInSpaceTime is base of interpolation calculated
	 * @return vector with interpolated Riemann invariants
	 */
	PdeVector interpolateValuesAround(const int s, const Mesh& mesh,
			const RealD direction, const Iterator& it, const PdeVector& dx,
			const bool canInterpolateInSpaceTime) {
		outerInvariants.clear();
		PdeVector ans = PdeVector::Zeros();
		
		for (int k = 0; k < PdeVector::M; k++) {
			
			if (dx(k) == 0) {
			// special for exact hit
				ans(k) = mesh.pde(it)(k);
				continue;
			}
			
			// point to interpolate respectively to point by given iterator
			RealD shift = direction * dx(k);
			Cell t = mesh.findCellCrossedByTheRay(it, shift);
			RiemannInvariant u = 0;
			
			if (t.n == t.N) {
			// characteristic hits into body
			// second order interpolate inner value in triangle on current time layer
				u = interpolateInSpace(mesh, mesh.coordsD(it) + shift, t, k);
				
			} else if (t.n == 0) {
			// outer characteristic from border/contact node
				outerInvariants.push_back(k);
				
			} else if (t.n == t.N - 1) {
			// characteristic hits out of body going throughout border face
				if (canInterpolateInSpaceTime) {
					u = interpolateInSpaceTime(s, mesh, it, shift, t, k);
				} else {
					outerInvariants.push_back(k);
				}
				
			} else if (t.n == t.N - 2) {
			// exact hit to border edge(point)
				if (canInterpolateInSpaceTime) {
					u = interpolateInSpaceTime1D(s, mesh, it, shift, t, k);
				} else {
					outerInvariants.push_back(k);
				}
				
			}
			
			ans(k) = u;
		}
		
		return ans;
	}
	
	
	/** Interpolate invariant from space on current time layer (2D case) */
	RiemannInvariant interpolateInSpace(
			const Mesh& mesh, const Real2& query, const Cell& c, const int k) const {
		RiemannInvariantGradient g[CELL_POINTS_NUMBER];
		for (int i = 0; i < CELL_POINTS_NUMBER; i++) {
			g[i] = {gradients[mesh.getIndex(c(i))](0)(k),
			        gradients[mesh.getIndex(c(i))](1)(k)};
		}
		return TriangleInterpolator<RiemannInvariant>::hybridInterpolate(
				mesh.coordsD(c(0)), mesh.pde(c(0))(k), g[0],
				mesh.coordsD(c(1)), mesh.pde(c(1))(k), g[1],
				mesh.coordsD(c(2)), mesh.pde(c(2))(k), g[2],
				query);
	}
	
	
	/** Interpolate invariant from space on current time layer (3D case) */
	RiemannInvariant interpolateInSpace(
			const Mesh& mesh, const Real3& query, const Cell& c, const int k) const {
		RiemannInvariantGradient g[CELL_POINTS_NUMBER];
		for (int i = 0; i < CELL_POINTS_NUMBER; i++) {
			g[i] = {gradients[mesh.getIndex(c(i))](0)(k),
			        gradients[mesh.getIndex(c(i))](1)(k),
			        gradients[mesh.getIndex(c(i))](2)(k)};
		}
		return TetrahedronInterpolator<RiemannInvariant>::hybridInterpolate(
				mesh.coordsD(c(0)), mesh.pde(c(0))(k), g[0],
				mesh.coordsD(c(1)), mesh.pde(c(1))(k), g[1],
				mesh.coordsD(c(2)), mesh.pde(c(2))(k), g[2],
				mesh.coordsD(c(3)), mesh.pde(c(3))(k), g[3],
				query);
	}
	
	
	/** 
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border in some point. It's possible either for border and inner nodes.
	 * @note border nodes must be already calculated
	 */
	RiemannInvariant interpolateInSpaceTime(const int s, const Mesh& mesh, 
			const Iterator& it, const Real2& shift, const Cell& borderEdge, const int k) const {
		/// 2D case
		/// first order interpolate in triangle formed by border points from
		/// current and next time layers (triangle in space-time)

		Real2 r1 = mesh.coordsD(borderEdge(0));
		Real2 r2 = mesh.coordsD(borderEdge(1));
		Real2 r0 = mesh.coordsD(it);
		Real2 rc = linal::linesIntersection(r1, r2, r0, r0 + shift);
				//< coordinate of border-characteristic intersection
		
		return TriangleInterpolator<RiemannInvariant>::interpolateInOwner(
				// current time layer
				{0, 0}, mesh.pde(borderEdge(0))(k),
				{1, 0}, mesh.pde(borderEdge(1))(k),
				// next time layer
				{0, 1}, mesh.pdeNew(s, borderEdge(0))(k),
				{1, 1}, mesh.pdeNew(s, borderEdge(1))(k),
				// query in space-time
				{    linal::length(rc - r1) / linal::length(r2 - r1),
				 1 - linal::length(rc - r0) / linal::length(shift)});
	}
	
	
	/** 
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border in some point. It's possible either for border and inner nodes.
	 * @note border nodes must be already calculated
	 */
	RiemannInvariant interpolateInSpaceTime(const int s, const Mesh& mesh, 
			const Iterator& it, const Real3& shift, const Cell& borderFace, const int k) const {
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
		return TetrahedronInterpolator<RiemannInvariant>::interpolateInOwner(
				// current time layer
				{0, 0, 0}, mesh.pde(borderFace(0))(k),
				{1, 0, 0}, mesh.pde(borderFace(1))(k),
				{0, 1, 0}, mesh.pde(borderFace(2))(k),
				// next time layer
				{0, 0, 1}, mesh.pdeNew(s, borderFace(0))(k),
				{1, 0, 1}, mesh.pdeNew(s, borderFace(1))(k),
				{0, 1, 1}, mesh.pdeNew(s, borderFace(2))(k),
				// query in space-time
				{w(0), w(1), 1 - linal::length(rc - r0) / linal::length(shift)});
	}
	
	
	/** 
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border exactly through the border vertex. It's possible for 2D case.
	 * It's possible either for border and inner nodes.
	 * @note border nodes must be already calculated
	 */
	RiemannInvariant interpolateInSpaceTime1D(const int s, const Mesh& mesh, 
			const Iterator& it, const Real2& shift, const Cell& borderVertex, const int k) const {
		/// 2D case
		/// first order interpolate in the line formed by crossed point
		/// at current and next time layers (line in space-time)

		auto bv = borderVertex(0);
		Real2 rv = mesh.coordsD(bv);
		Real2 r0 = mesh.coordsD(it);
		real w = linal::length(rv - r0) / linal::length(shift);
		return mesh.pde(bv)(k) * w + mesh.pdeNew(s, it)(k) * (1 - w);
	}
	
	
	/** 
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border exactly through the border vertex. It's possible for 2D case.
	 * It's possible either for border and inner nodes.
	 * @note border nodes must be already calculated
	 */
	RiemannInvariant interpolateInSpaceTime1D(const int /*s*/, const Mesh& /*mesh*/, 
			const Iterator& /*it*/, const Real3& /*shift*/, const Cell& /*borderEdge*/, const int /*k*/) const {
		THROW_UNSUPPORTED("TODO");
	}
	
	
	/// List of outer Riemann invariants after node calculation.
	/// Invariants are specified by their indices in matrix L.
	WaveIndices outerInvariants;
	
	/// The additional storage of PDE-vectors for the opportunity
	/// to save some time layers temporary during stage calculation
	std::vector<PdeVariables> savedPdeTimeLayer;
	
	/// The storage of gradients of mesh pde values.
	std::vector<PdeGradient> gradients;
	/// The storage of hessians of mesh pde values. (Unused now)
	std::vector<PdeHessian> hessians;
	
	USE_AND_INIT_LOGGER("gcm.simplex.GridCharacteristicMethod")
};


} // namespace simplex
} // namespace gcm


#endif // LIBGCM_SIMPLEX_GRIDCHARACTERISTICMETHODINRIEMANNINVARIANTS_HPP
