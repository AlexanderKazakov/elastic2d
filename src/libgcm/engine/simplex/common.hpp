#ifndef LIBGCM_SIMPLEX_COMMON_HPP
#define LIBGCM_SIMPLEX_COMMON_HPP

#include <libgcm/grid/AbstractGrid.hpp>
#include <libgcm/util/math/GridCharacteristicMethod.hpp>
#include <libgcm/util/math/Differentiation.hpp>
#include <libgcm/util/math/interpolation/interpolation.hpp>


namespace gcm {
namespace simplex {

/**
 * Interface for different types of gcm-methods
 */
class GridCharacteristicMethodBase {
public:
	/** Actions necessary before the stage */
	virtual void beforeStage(
			const int nextPdeLayerIndex,
			const int s, AbstractGrid& mesh_) = 0;
	
	/**
	 * Do inner part of stage of splitting grid-characteristic method
	 * on contact and border nodes
	 */
	virtual void contactAndBorderStage(
			const int nextPdeLayerIndex,
			const int s, const real timeStep, AbstractGrid& mesh_) = 0;
	
	/**
	 * Do stage of splitting grid-characteristic method on inner nodes
	 */
	virtual void innerStage(
			const int nextPdeLayerIndex,
			const int s, const real timeStep, AbstractGrid& mesh_) = 0;
	
	/** Actions necessary after the stage */
	virtual void afterStage(
			const int nextPdeLayerIndex,
			const int s, AbstractGrid& mesh_) = 0;
	
	
protected:
	/** Points where characteristics from next time layer cross current time layer */
	template<typename Mesh>
	static inline
	typename Mesh::PdeVector crossingPoints(
			const typename Mesh::Iterator& it, const int s,
			const real timeStep, const Mesh& mesh) {
		return -timeStep * linal::diag(mesh.matrices(it)->m[s].L);
	}
	
	
	/**
	 * Handle the case when characteristic goes inside the body and then cross 
	 * the border in some point (2D case):
	 * first order interpolate in triangle formed by border points from
	 * current and next time layers (triangle in space-time)
	 * @note border nodes must be already calculated
	 * @param r0 coordinates of node in calculation
	 * @param r_i coordinates of two border nodes
	 * @param v_i_curr values in two border nodes on current time layer
	 * @param v_i_next values in two border nodes on next time layer
	 * @return space-time interpolated value at the point of
	 * border and characteristic intersection
	 */
	template<typename Value>
	static inline
	Value interpolateInSpaceTime(Real2 shift, Real2 r0,
			Real2 r1, Real2 r2,
			Value v1curr, Value v2curr,
			Value v1next, Value v2next) {
		// coordinate of border with characteristic intersection
		Real2 rc = linal::linesIntersection(r1, r2, r0, r0 + shift);
		return TriangleInterpolator<Value>::interpolateInOwner(
				// current time layer
				{0, 0}, v1curr,
				{1, 0}, v2curr,
				// next time layer
				{0, 1}, v1next,
				{1, 1}, v2next,
				// query in space-time
				{    linal::length(rc - r1) / linal::length(r2 - r1),
				 1 - linal::length(rc - r0) / linal::length(shift)});
	}
	
	
	/**
	 * Handle the case when characteristic goes inside the body and then cross 
	 * the border in some point (3D case):
	 * first order interpolate in tetrahedron formed by border points from
	 * current and next time layers (tetrahedron in space-time)
	 * @note border nodes must be already calculated
	 * @param r0 coordinates of node in calculation
	 * @param r_i coordinates of three border nodes
	 * @param v_i_curr values in three border nodes on current time layer
	 * @param v_i_next values in three border nodes on next time layer
	 * @return space-time interpolated value at the point of
	 * border and characteristic intersection
	 */
	template<typename Value>
	static inline
	Value interpolateInSpaceTime(Real3 shift, Real3 r0,
			Real3 r1, Real3 r2, Real3 r3,
			Value v1curr, Value v2curr, Value v3curr,
			Value v1next, Value v2next, Value v3next) {
		// coordinate of border-characteristic intersection
		Real3 rc = linal::lineWithFlatIntersection(r1, r2, r3, r0, r0 + shift);
		// find local coordinates of (rc - r1) in basis {e1, e2}
		Real3 e1 = r2 - r1, e2 = r3 - r1, a = rc - r1;
		linal::Matrix<3, 2> A = {
				e1(0), e2(0),
				e1(1), e2(1),
				e1(2), e2(2),
		};
		Real2 w = linal::linearLeastSquares(A, a);
		return TetrahedronInterpolator<Value>::interpolateInOwner(
				// current time layer
				{0, 0, 0}, v1curr,
				{1, 0, 0}, v2curr,
				{0, 1, 0}, v3curr,
				// next time layer
				{0, 0, 1}, v1next,
				{1, 0, 1}, v2next,
				{0, 1, 1}, v3next,
				// query in space-time
				{w(0), w(1), 1 - linal::length(rc - r0) / linal::length(shift)});
	}
	
	
	/**
	 * Handle the case when characteristic goes inside the body and then cross 
	 * border exactly through the border vertex (2D)
	 * @note border nodes must be already calculated
	 */
	template<typename Value>
	static inline
	Value interpolateInSpaceTime1D(Real2 shift, Real2 r0,
			Real2 r1, Value v1curr, Value v1next) {
		/// first order interpolate in the line formed by crossed point
		/// at current and next time layers (segment in space-time)
		real w = linal::length(r1 - r0) / linal::length(shift);
		return v1curr * w + v1next * (1 - w);
	}
};


/**
 * Get specified columns from gcmMatrices at specified stage.
 * Useful for border/contact correctors to get outer invariants.
 */
template<typename Model>
inline
typename Model::OuterMatrix getColumnsFromGcmMatrices(
		const int stage,
		const typename Model::WaveIndices columnsIndices,
		std::shared_ptr<const typename Model::GCM_MATRICES> m) {
	typename Model::OuterMatrix ans = Model::OuterMatrix::Zeros();
	int counter = 0;
	for (const int i : columnsIndices) {
		ans.setColumn(counter++, (*m)(stage).U1.getColumn(i));
	}
	return ans;
}


} // namespace simplex
} // namespace gcm


#endif // LIBGCM_SIMPLEX_COMMON_HPP
