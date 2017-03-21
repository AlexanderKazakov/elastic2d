#ifndef LIBGCM_SIMPLEX_COMMON_HPP
#define LIBGCM_SIMPLEX_COMMON_HPP

#include <libgcm/util/infrastructure/infrastructure.hpp>
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
