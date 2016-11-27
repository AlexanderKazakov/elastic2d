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
	virtual void beforeStage(const AbstractGrid& mesh_) = 0;
	virtual void contactStage(
			const int s, const real timeStep, AbstractGrid& mesh_) = 0;
	virtual void stage(
			const int s, const real timeStep, AbstractGrid& mesh_) = 0;
	virtual void returnBackBadOuterCases(
			const int s, AbstractGrid& mesh_) const = 0;
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
	static const int DIMENSIONALITY = Mesh::DIMENSIONALITY;
	
	virtual void beforeStage(const AbstractGrid& mesh_) override {
		/// calculate spatial derivatives of all mesh pde values ones before stage
		/// in order to use them multiple times while stage calculation 
		const Mesh& mesh = dynamic_cast<const Mesh&>(mesh_);
		DIFFERENTIATION::estimateGradient(mesh, gradients);
	}
	
	
	/**
	 * Do grid-characteristic stage of splitting method on contact and border nodes
	 */
	virtual void contactStage(
			const int s, const real timeStep, AbstractGrid& mesh_) override {
		Mesh& mesh = dynamic_cast<Mesh&>(mesh_);
		RealD direction = mesh.getCalculationBasis().getColumn(s);
		outerCasesToReturnBack.clear();
		
		/// calculate inner waves of contact nodes
//		#pragma omp parallel for
		for (auto contactIter = mesh.contactBegin(); 
		          contactIter < mesh.contactEnd(); ++contactIter) {
			mesh._pdeNew(s, *contactIter) = localGcmStep(
					mesh.matrices(*contactIter)->m[s].U1,
					mesh.matrices(*contactIter)->m[s].U,
					interpolateValuesAround(s, mesh, direction, *contactIter,
							crossingPoints(*contactIter, s, timeStep, mesh), false));
			checkOuterCases(s, *contactIter, mesh);
		}
		
		/// calculate inner waves of border nodes
//		#pragma omp parallel for
		for (auto borderIter = mesh.borderBegin(); 
		          borderIter < mesh.borderEnd(); ++borderIter) {
			mesh._pdeNew(s, *borderIter) = localGcmStep(
					mesh.matrices(*borderIter)->m[s].U1,
					mesh.matrices(*borderIter)->m[s].U,
					interpolateValuesAround(s, mesh, direction, *borderIter,
							crossingPoints(*borderIter, s, timeStep, mesh), false));
			checkOuterCases(s, *borderIter, mesh);
		}
	}
	
	
	/**
	 * Do grid-characteristic stage of splitting method on inner nodes
	 * @note contact and border nodes must be already calculated
	 * @param s number of stage (GcmMatrix number)
	 * @param timeStep time step
	 * @param mesh mesh to perform calculation
	 */
	virtual void stage(
			const int s, const real timeStep, AbstractGrid& mesh_) override {
		Mesh& mesh = dynamic_cast<Mesh&>(mesh_);
		RealD direction = mesh.getCalculationBasis().getColumn(s);
		
		/// calculate inner nodes
//		#pragma omp parallel for
		for (auto innerIter = mesh.innerBegin(); 
		          innerIter < mesh.innerEnd(); ++innerIter) {
			mesh._pdeNew(s, *innerIter) = localGcmStep(
					mesh.matrices(*innerIter)->m[s].U1,
					mesh.matrices(*innerIter)->m[s].U,
					interpolateValuesAround(s, mesh, direction, *innerIter,
							crossingPoints(*innerIter, s, timeStep, mesh), true));
			assert_eq(outerInvariants.size(), 0);
		}
	}
	
	
	/**
	 * Rewrite back pdeNew values, corrected by border and contact correctors,
	 * because for double-outer cases such correction is invalid.
	 */
	virtual void returnBackBadOuterCases(
			const int s, AbstractGrid& mesh_) const override {
		Mesh& mesh = dynamic_cast<Mesh&>(mesh_);
		for (std::pair<Iterator, PdeVector> p : outerCasesToReturnBack) {
			mesh._pdeNew(s, p.first) = p.second;
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
	 * @param s stage
	 * @param mesh mesh to perform interpolation on
	 * @param direction direction of line to find values along
	 * @param it index-iterator of node
	 * @param dx Vector of distances from reference node on which
	 * values should be interpolated
	 * @param canInterpolateInSpaceTime is base of interpolation calculated
	 * @return Matrix with interpolated nodal values in columns
	 */
	Matrix interpolateValuesAround(const int s, 
	                               const Mesh& mesh, const RealD direction,
	                               const Iterator& it, const PdeVector& dx,
	                               const bool canInterpolateInSpaceTime) {
		outerInvariants.clear();
		Matrix ans = Matrix::Zeros();
		
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
				
			} else if (t.n == 0) {
			// outer characteristic from border/contact node
				outerInvariants.push_back(k);
				
			} else if (t.n == t.N - 1 && canInterpolateInSpaceTime) {
			// characteristic hits out of body going throughout border face
				u = interpolateInSpaceTime(s, mesh, it, shift, t);
				
			} else if (t.n == t.N - 2 && canInterpolateInSpaceTime) {
			// exact hit to border edge(point)
				u = interpolateInSpaceTime1D(s, mesh, it, shift, t);
				
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
	PdeVector interpolateInSpaceTime(const int s, const Mesh& mesh, 
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
				{0, 1}, mesh.pdeNew(s, borderEdge(0)),
				{1, 1}, mesh.pdeNew(s, borderEdge(1)),
				// query in space-time
				{    linal::length(rc - r1) / linal::length(r2 - r1),
				 1 - linal::length(rc - r0) / linal::length(shift)});
	}
	
	
	/** 
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border in some point. It's possible either for border and inner nodes.
	 * @note border nodes must be already calculated
	 */
	PdeVector interpolateInSpaceTime(const int s, const Mesh& mesh, 
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
				{0, 0, 1}, mesh.pdeNew(s, borderFace(0)),
				{1, 0, 1}, mesh.pdeNew(s, borderFace(1)),
				{0, 1, 1}, mesh.pdeNew(s, borderFace(2)),
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
	PdeVector interpolateInSpaceTime1D(const int s, const Mesh& mesh, 
			const Iterator& it, const Real2& shift, const Cell& borderVertex) const {
		/// 2D case
		/// first order interpolate in the line formed by crossed point
		/// at current and next time layers (line in space-time)

		auto bv = borderVertex(0);
		Real2 rv = mesh.coordsD(bv);
		Real2 r0 = mesh.coordsD(it);
		real w = linal::length(rv - r0) / linal::length(shift);
		return mesh.pde(bv) * w + mesh.pdeNew(s, it) * (1 - w);
	}
	
	
	/** 
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border exactly through the border vertex. It's possible for 2D case.
	 * It's possible either for border and inner nodes.
	 * @note border nodes must be already calculated
	 */
	PdeVector interpolateInSpaceTime1D(const int /*s*/, const Mesh& /*mesh*/, 
			const Iterator& /*it*/, const Real3& /*shift*/, const Cell& /*borderEdge*/) const {
		THROW_UNSUPPORTED("TODO");
		return PdeVector::Zeros(); // FIXME
	}
	
	
	/**
	 * Check and remember cases with "illegal" outer characteristics
	 */
	void checkOuterCases(const int s, const Iterator it, const Mesh& mesh) {
		
		if (outerInvariants == Model::LEFT_INVARIANTS ||
			outerInvariants == Model::RIGHT_INVARIANTS) {
		/// The normal case for border corrector
			return;
		}
		
		if (outerInvariants.size() == 0) {
		/// No outers. Nothing to do for border corrector
			outerCasesToReturnBack.push_back({it, mesh.pdeNew(s, it)});
			return;
		}
		
		if (outerInvariants.size() == 2 * OUTER_NUMBER) {
		/// All outers. No space -- no waves. Set to old value
			outerCasesToReturnBack.push_back({it, mesh.pde(it)});
			return;
		}
		
		/// Some geometrical inexactness. Most likely that the calculation
		/// direction is almost parallel to the border. 
		/// So disable border correction
		/// TODO - make special for such cases more precise search?
//		std::cout << outerInvariants.size() << " outer characteristics: ";
//		for (int k : outerInvariants) { std::cout << k << " "; }
//		std::cout << " at:" << mesh.coordsD(it);
		outerCasesToReturnBack.push_back({it, mesh.pdeNew(s, it)});
	}
	
	
	/// List of outer Riemann invariants after node calculation.
	/// Invariants are specified by their indices in matrix L.
	std::vector<int> outerInvariants;
	
	/// List of nodes which have "illegal" outer characteristics on the current
	/// calculation direction. This is possible for borders and contacts only.
	std::vector<std::pair<Iterator, PdeVector>> outerCasesToReturnBack;
	
	/// The storage of gradients of mesh pde values.
	std::vector<PdeGradient> gradients;
	/// The storage of hessians of mesh pde values. (Unused now)
	std::vector<PdeHessian> hessians;
	
	USE_AND_INIT_LOGGER("gcm.simplex.GridCharacteristicMethod")
	
};


} // namespace simplex
} // namespace gcm


#endif // LIBGCM_SIMPLEX_GRIDCHARACTERISTICMETHOD_HPP
