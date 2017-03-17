#ifndef LIBGCM_SIMPLEX_COMMON_HPP
#define LIBGCM_SIMPLEX_COMMON_HPP

namespace gcm {
namespace simplex {

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
