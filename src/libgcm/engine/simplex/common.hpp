#ifndef LIBGCM_SIMPLEX_COMMON_HPP
#define LIBGCM_SIMPLEX_COMMON_HPP

namespace gcm {
namespace simplex {

template<typename Model>
inline
typename Model::OuterMatrix getOuterMatrixFromGcmMatricesInLocalBasis(
		std::shared_ptr<const typename Model::GCM_MATRICES> m) {
	typename Model::OuterMatrix ans = Model::OuterMatrix::Zeros();
	int counter = 0;
	for (const int i : Model::RIGHT_INVARIANTS) {
		ans.setColumn(counter++, (*m)(0).U1.getColumn(i));
	}
	return ans;
}

} // namespace simplex
} // namespace gcm


#endif // LIBGCM_SIMPLEX_COMMON_HPP
