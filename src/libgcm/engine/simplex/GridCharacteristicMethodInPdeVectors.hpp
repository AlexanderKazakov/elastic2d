#ifndef LIBGCM_SIMPLEX_GRIDCHARACTERISTICMETHODINPDEVECTORS_HPP
#define LIBGCM_SIMPLEX_GRIDCHARACTERISTICMETHODINPDEVECTORS_HPP

#include <libgcm/engine/simplex/common.hpp>

namespace gcm {
namespace simplex {

/**
 * Grid-characteristic method for meshes based on SimplexGrid.
 * The approach is to calculate and advect along characteristics
 * PDE-vectors not scalar Riemann-invariants
 * @see GridCharacteristicMethodInRiemannInvariants -- an opposite approach
 */
template<typename Mesh>
class GridCharacteristicMethodInPdeVectors :
		public GridCharacteristicMethodBase {
public:
	typedef GridCharacteristicMethodBase                       Base;
	
	typedef typename Mesh::Matrix                              Matrix;
	typedef typename Mesh::PdeVector                           PdeVector;
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
	
	virtual void beforeStage(
			const int /*nextPdeLayerIndex*/,
			const int /*s*/, AbstractGrid& mesh_) override {
		/// calculate spatial derivatives of all mesh pde values ones before stage
		/// in order to use them multiple times while stage calculation 
		const Mesh& mesh = dynamic_cast<const Mesh&>(mesh_);
		DIFFERENTIATION::estimateGradient(mesh, gradients);
	}
	
	
	virtual void contactAndBorderStage(
			const int nextPdeLayerIndex,
			const int s, const real timeStep, AbstractGrid& mesh_) override {
		Mesh& mesh = dynamic_cast<Mesh&>(mesh_);
		
		/// calculate inner waves of contact nodes
		for (auto iter = mesh.contactBegin(); iter < mesh.contactEnd(); ++iter) {
			const GcmMatrices& gcmMatrices = *mesh.matrices(*iter);
			const RealD direction = gcmMatrices.basis.getColumn(s);
			mesh._pdeNew(nextPdeLayerIndex, *iter) = localGcmStep(
				gcmMatrices(s).U1, gcmMatrices(s).U,
				interpolateValuesAround(nextPdeLayerIndex, s, mesh, direction, *iter,
					Base::crossingPoints(*iter, s, timeStep, mesh), false));
			mesh._waveIndices(*iter) = outerInvariants;
		}
		
		/// calculate inner waves of border nodes
		for (auto iter = mesh.borderBegin(); iter < mesh.borderEnd(); ++iter) {
			const GcmMatrices& gcmMatrices = *mesh.matrices(*iter);
			const RealD direction = gcmMatrices.basis.getColumn(s);
			mesh._pdeNew(nextPdeLayerIndex, *iter) = localGcmStep(
				gcmMatrices(s).U1, gcmMatrices(s).U,
				interpolateValuesAround(nextPdeLayerIndex, s, mesh, direction, *iter,
					Base::crossingPoints(*iter, s, timeStep, mesh), false));
			mesh._waveIndices(*iter) = outerInvariants;
		}
	}
	
	
	virtual void innerStage(
			const int nextPdeLayerIndex,
			const int s, const real timeStep, AbstractGrid& mesh_) override {
		Mesh& mesh = dynamic_cast<Mesh&>(mesh_);
		const RealD direction = mesh.getInnerCalculationBasis().getColumn(s);
		
		/// calculate inner nodes
		for (auto innerIter = mesh.innerBegin(); 
		          innerIter < mesh.innerEnd(); ++innerIter) {
			mesh._pdeNew(nextPdeLayerIndex, *innerIter) = localGcmStep(
				mesh.matrices(*innerIter)->m[s].U1,
				mesh.matrices(*innerIter)->m[s].U,
				interpolateValuesAround(nextPdeLayerIndex, s, mesh, direction, *innerIter,
					Base::crossingPoints(*innerIter, s, timeStep, mesh), true));
			assert_eq(outerInvariants.size(), 0);
		}
	}
	
	
	virtual void afterStage(
			const int /*nextPdeLayerIndex*/,
			const int /*s*/, AbstractGrid& /*mesh_*/) override { }
	
	
private:
	/**
	 * Interpolate nodal values in specified points.
	 * Interpolated value for k-th point in vector dx are
	 * stored in k-th column of returned Matrix.
	 * If specified point appears to be out of body
	 * AND it is really border case, matrix column is set to zeros
	 * and outerInvariants is added with the index.
	 * @param nextPdeLayerIndex not always equal to stage!
	 * @param s stage
	 * @param mesh mesh to perform interpolation on
	 * @param direction direction of line to find values along
	 * @param it index-iterator of node
	 * @param dx Vector of distances from reference node on which
	 * values should be interpolated
	 * @param canInterpolateInSpaceTime is base of interpolation calculated
	 * @return Matrix with interpolated nodal values in columns
	 */
	Matrix interpolateValuesAround(const int nextPdeLayerIndex, const int s, 
	                               const Mesh& mesh, const RealD direction,
	                               const Iterator& it, const PdeVector& dx,
	                               const bool canInterpolateInSpaceTime) {
		outerInvariants.clear();
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
			// characteristic hits into body
			// second order interpolate inner value in triangle on current time layer
				u = interpolateInSpace(mesh, mesh.coordsD(it) + shift, t);
				
			} else if (t.n == 0) {
			// outer characteristic from border/contact node
				outerInvariants.push_back(k);
				
			} else if (t.n == t.N - 1) {
			// characteristic hits out of body going throughout border face
				if (canInterpolateInSpaceTime) {
					u = interpolateInSpaceTime(nextPdeLayerIndex, s, mesh, it, shift, t);
				} else {
					outerInvariants.push_back(k);
				}
				
			} else if (t.n == t.N - 2) {
			// exact hit to border edge(point)
				if (canInterpolateInSpaceTime) {
					u = interpolateInSpaceTime1D(nextPdeLayerIndex, s, mesh, it, shift, t);
				} else {
					outerInvariants.push_back(k);
				}
				
			}
			
			ans.setColumn(k, u);
		}
		
		return ans;
	}
	
	
	/** Interpolate PdeVector from space on current time layer (2D case) */
	inline
	PdeVector interpolateInSpace(const Mesh& mesh, const Real2& query, const Cell& c) const {
		return TriangleInterpolator<PdeVector>::hybridInterpolate(
				mesh.coordsD(c(0)), mesh.pde(c(0)), gradients[mesh.getIndex(c(0))],
				mesh.coordsD(c(1)), mesh.pde(c(1)), gradients[mesh.getIndex(c(1))],
				mesh.coordsD(c(2)), mesh.pde(c(2)), gradients[mesh.getIndex(c(2))],
				query);
	}
	
	
	/** Interpolate PdeVector from space on current time layer (3D case) */
	inline
	PdeVector interpolateInSpace(const Mesh& mesh, const Real3& query, const Cell& c) const {
		return TetrahedronInterpolator<PdeVector>::hybridInterpolate(
				mesh.coordsD(c(0)), mesh.pde(c(0)), gradients[mesh.getIndex(c(0))],
				mesh.coordsD(c(1)), mesh.pde(c(1)), gradients[mesh.getIndex(c(1))],
				mesh.coordsD(c(2)), mesh.pde(c(2)), gradients[mesh.getIndex(c(2))],
				mesh.coordsD(c(3)), mesh.pde(c(3)), gradients[mesh.getIndex(c(3))],
				query);
	}
	
	
	/** 
	 * Handle the case when characteristic goes inside the body
	 * and then cross the border in some point (2D)
	 * @note border nodes must be already calculated
	 */
	static inline
	PdeVector interpolateInSpaceTime(
			const int nextPdeLayerIndex, const int /*s*/, const Mesh& mesh, 
			const Iterator& it, const Real2& shift, const Cell& borderEdge) {
		return Base::interpolateInSpaceTime(shift, mesh.coordsD(it),
				mesh.coordsD(borderEdge(0)), mesh.coordsD(borderEdge(1)),
				mesh.pde(borderEdge(0)), mesh.pde(borderEdge(1)),
				mesh.pdeNew(nextPdeLayerIndex, borderEdge(0)),
				mesh.pdeNew(nextPdeLayerIndex, borderEdge(1)));
	}
	
	
	/** 
	 * Handle the case when characteristic goes inside the body
	 * and then cross border in some point (3D)
	 * @note border nodes must be already calculated
	 */
	static inline
	PdeVector interpolateInSpaceTime(
			const int nextPdeLayerIndex, const int /*s*/, const Mesh& mesh, 
			const Iterator& it, const Real3& shift, const Cell& borderFace) {
		return Base::interpolateInSpaceTime(shift, mesh.coordsD(it),
				mesh.coordsD(borderFace(0)), mesh.coordsD(borderFace(1)),
				mesh.coordsD(borderFace(2)),
				mesh.pde(borderFace(0)), mesh.pde(borderFace(1)),
				mesh.pde(borderFace(2)),
				mesh.pdeNew(nextPdeLayerIndex, borderFace(0)),
				mesh.pdeNew(nextPdeLayerIndex, borderFace(1)),
				mesh.pdeNew(nextPdeLayerIndex, borderFace(2)));
	}
	
	
	/**
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border exactly through the border vertex. (in 2D case)
	 * @note border nodes must be already calculated
	 */
	static inline
	PdeVector interpolateInSpaceTime1D(
			const int nextPdeLayerIndex, const int /*s*/, const Mesh& mesh, 
			const Iterator& it, const Real2& shift, const Cell& borderVertex) {
		return Base::interpolateInSpaceTime1D(shift, mesh.coordsD(it),
				mesh.coordsD(borderVertex(0)),
				mesh.pde(borderVertex(0)),
				mesh.pdeNew(nextPdeLayerIndex, borderVertex(0)));
	}
	
	
	/**
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border exactly through the border edge. (in 3D)
	 * @note border nodes must be already calculated
	 */
	static inline
	PdeVector interpolateInSpaceTime1D(
			const int /*nextPdeLayerIndex*/, const int /*s*/, const Mesh& /*mesh*/, 
			const Iterator& /*it*/, const Real3& /*shift*/, const Cell& /*borderVertex*/) {
		THROW_UNSUPPORTED("This was not occured ever before");
	}
	
	
	/// List of outer Riemann invariants after node calculation.
	/// Invariants are specified by their indices in matrix L.
	WaveIndices outerInvariants;
	
	/// The storage of gradients of mesh pde values.
	std::vector<PdeGradient> gradients;
	/// The storage of hessians of mesh pde values. (Unused now)
	std::vector<PdeHessian> hessians;
	
	USE_AND_INIT_LOGGER("gcm.simplex.GridCharacteristicMethodInPdeVectors")
};


} // namespace simplex
} // namespace gcm


#endif // LIBGCM_SIMPLEX_GRIDCHARACTERISTICMETHODINPDEVECTORS_HPP
