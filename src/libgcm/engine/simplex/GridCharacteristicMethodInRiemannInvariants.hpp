#ifndef LIBGCM_SIMPLEX_GRIDCHARACTERISTICMETHODINRIEMANNINVARIANTS_HPP
#define LIBGCM_SIMPLEX_GRIDCHARACTERISTICMETHODINRIEMANNINVARIANTS_HPP

#include <libgcm/engine/simplex/common.hpp>

namespace gcm {
namespace simplex {

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
	typedef GridCharacteristicMethodBase                       Base;
	
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
	
	typedef real                                  RiemannInvariant;
	typedef linal::VECTOR<
			DIMENSIONALITY, RiemannInvariant>     RiemannInvariantGradient;
	
	
	virtual void beforeStage(
			const int /*nextPdeLayerIndex*/,
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
	
	
	virtual void contactAndBorderStage(
			const int nextPdeLayerIndex, const int s,
			const real timeStep, AbstractGrid& mesh_) override {
		Mesh& mesh = dynamic_cast<Mesh&>(mesh_);
		
		auto calculate = [&](const Iterator iter) {
			const GcmMatrices& gcmMatrices = *mesh.matrices(iter);
			const RealD direction = gcmMatrices.basis.getColumn(s);
			mesh._pdeNew(nextPdeLayerIndex, iter) = interpolateValuesAround(
					nextPdeLayerIndex,
					s, mesh, direction, iter,
					Base::crossingPoints(iter, s, timeStep, mesh), false);
			
			if (outerInvariants != Model::RIGHT_INVARIANTS &&
					outerInvariants != Model::LEFT_INVARIANTS &&
					outerInvariants.size() != 2 * Model::OUTER_NUMBER &&
					!outerInvariants.empty()) {
				if (!Utils::intersection(outerInvariants, Model::RIGHT_INVARIANTS).empty()) {
					outerInvariants = Utils::summ(outerInvariants, Model::RIGHT_INVARIANTS);
				}
				if (!Utils::intersection(outerInvariants, Model::LEFT_INVARIANTS).empty()) {
					outerInvariants = Utils::summ(outerInvariants, Model::LEFT_INVARIANTS);
				}
//				for (int i : outerInvariants) {
//					mesh._pdeNew(nextPdeLayerIndex, iter)(i) = 0;
//				}
			}
			
			mesh._waveIndices(iter) = outerInvariants;
		};
		
		for (auto iter = mesh.contactBegin(); iter < mesh.contactEnd(); ++iter) {
			calculate(*iter);
		}
		for (auto iter = mesh.borderBegin(); iter < mesh.borderEnd(); ++iter) {
			calculate(*iter);
		}
	}
	
	
	virtual void innerStage(
			const int nextPdeLayerIndex, const int s,
			const real timeStep, AbstractGrid& mesh_) override {
		Mesh& mesh = dynamic_cast<Mesh&>(mesh_);
		const RealD direction = mesh.getInnerCalculationBasis().getColumn(s);
		
		/// calculate inner nodes
		for (auto iter = mesh.innerBegin(); iter < mesh.innerEnd(); ++iter) {
			mesh._pdeNew(nextPdeLayerIndex, *iter) = interpolateValuesAround(
					nextPdeLayerIndex,
					s, mesh, direction, *iter,
					Base::crossingPoints(*iter, s, timeStep, mesh), true);
			assert_eq(outerInvariants.size(), 0);
		}
	}
	
	
	virtual void afterStage(
			const int nextPdeLayerIndex,
			const int s, AbstractGrid& mesh_) override {
		Mesh& mesh = dynamic_cast<Mesh&>(mesh_);
		/// set PDE-vectors on old time layer to its values saved before the stage
		mesh.swapPdeVariablesStorage(savedPdeTimeLayer);
		/// switch back from Riemann invariants to PDE-variables
		for (const Iterator it : mesh) {
			mesh._pdeNew(nextPdeLayerIndex, it) =
					(*mesh.matrices(it))(s).U1 * mesh.pdeNew(nextPdeLayerIndex, it);
		}
	}
	
	
private:
	/**
	 * Interpolate Riemann invariants in specified points.
	 * If specified point appears to be out of body
	 * AND it is really border case, invariant is set to zero
	 * and outerInvariants is added with the index.
	 * @param nextPdeLayerIndex not always equal to stage!
	 * @param s stage
	 * @param mesh mesh to perform interpolation on
	 * @param direction direction of line to find values along
	 * @param it index-iterator of node
	 * @param dx Vector of distances from reference node on which
	 * values should be interpolated
	 * @param canInterpolateInSpaceTime is base of interpolation calculated
	 * @return vector with interpolated Riemann invariants
	 */
	PdeVector interpolateValuesAround(
			const int nextPdeLayerIndex,
			const int s,const Mesh& mesh,
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
					u = interpolateInSpaceTime(s, mesh, it, shift, t, nextPdeLayerIndex, k);
				} else {
					outerInvariants.push_back(k);
				}
				
			} else if (t.n == t.N - 2) {
			// exact hit to border edge(point)
				if (canInterpolateInSpaceTime) {
					u = interpolateInSpaceTime1D(s, mesh, it, shift, t, nextPdeLayerIndex, k);
				} else {
					outerInvariants.push_back(k);
				}
				
			}
			
			ans(k) = u;
		}
		
		return ans;
	}
	
	
	/** Interpolate invariant from space on current time layer (2D case) */
	inline
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
	inline
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
	 * Handle the case when characteristic goes inside the body
	 * and then cross the border in some point (2D)
	 * @note border nodes must be already calculated
	 */
	static inline
	RiemannInvariant interpolateInSpaceTime(const int /*s*/, const Mesh& mesh, 
			const Iterator& it, const Real2& shift, const Cell& borderEdge,
			const int nextPdeLayerIndex, const int k) {
		return Base::interpolateInSpaceTime(shift, mesh.coordsD(it),
				mesh.coordsD(borderEdge(0)), mesh.coordsD(borderEdge(1)),
				mesh.pde(borderEdge(0))(k), mesh.pde(borderEdge(1))(k),
				mesh.pdeNew(nextPdeLayerIndex, borderEdge(0))(k),
				mesh.pdeNew(nextPdeLayerIndex, borderEdge(1))(k));
	}
	
	
	/**
	 * Handle the case when characteristic goes inside the body
	 * and then cross border in some point (3D)
	 * @note border nodes must be already calculated
	 */
	static inline
	RiemannInvariant interpolateInSpaceTime(const int /*s*/, const Mesh& mesh, 
			const Iterator& it, const Real3& shift, const Cell& borderFace,
			const int nextPdeLayerIndex, const int k) {
		return Base::interpolateInSpaceTime(shift, mesh.coordsD(it),
				mesh.coordsD(borderFace(0)), mesh.coordsD(borderFace(1)),
				mesh.coordsD(borderFace(2)),
				mesh.pde(borderFace(0))(k), mesh.pde(borderFace(1))(k),
				mesh.pde(borderFace(2))(k),
				mesh.pdeNew(nextPdeLayerIndex, borderFace(0))(k),
				mesh.pdeNew(nextPdeLayerIndex, borderFace(1))(k),
				mesh.pdeNew(nextPdeLayerIndex, borderFace(2))(k));
	}
	
	
	/**
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border exactly through the border vertex. (in 2D case)
	 * @note border nodes must be already calculated
	 */
	static inline
	RiemannInvariant interpolateInSpaceTime1D(const int /*s*/, const Mesh& mesh, 
			const Iterator& it, const Real2& shift, const Cell& borderVertex,
			const int nextPdeLayerIndex, const int k) {
		return Base::interpolateInSpaceTime1D(shift, mesh.coordsD(it),
				mesh.coordsD(borderVertex(0)),
				mesh.pde(borderVertex(0))(k),
				mesh.pdeNew(nextPdeLayerIndex, borderVertex(0))(k));
	}
	
	
	/**
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border exactly through the border edge. (in 3D)
	 * @note border nodes must be already calculated
	 */
	static inline
	RiemannInvariant interpolateInSpaceTime1D(const int /*s*/, const Mesh& /*mesh*/, 
			const Iterator& /*it*/, const Real3& /*shift*/,
			const Cell& /*borderVertex*/, const int /*nextPdeLayerIndex*/, const int /*k*/) {
		THROW_UNSUPPORTED("This did not occur ever before");
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
	
	USE_AND_INIT_LOGGER("gcm.simplex.GridCharacteristicMethodInRiemannInvariants")
};


} // namespace simplex
} // namespace gcm


#endif // LIBGCM_SIMPLEX_GRIDCHARACTERISTICMETHODINRIEMANNINVARIANTS_HPP
