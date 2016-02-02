#include <lib/rheology/gcm_matrices/GcmMatrices.hpp>
#include <lib/rheology/variables/variables.hpp>
#include <lib/rheology/materials/materials.hpp>

using namespace gcm;


template<>
const std::map<Waves::T, int/* number of column in U1 */> GcmMatrices<VelocitySigmaVariables<1>, IsotropicMaterial>::WAVE_COLUMNS = {
		{Waves::T::P_FORWARD,  0},
		{Waves::T::P_BACKWARD, 1}
};

template<>
const std::map<Waves::T, int/* number of column in U1 */> GcmMatrices<VelocitySigmaVariables<2>, IsotropicMaterial>::WAVE_COLUMNS = {
		{Waves::T::P_FORWARD,   1},
		{Waves::T::P_BACKWARD,  0},
		{Waves::T::S1_FORWARD,  3},
		{Waves::T::S1_BACKWARD, 2}
};

template<>
const std::map<Waves::T, int/* number of column in U1 */> GcmMatrices<VelocitySigmaVariables<3>, IsotropicMaterial>::WAVE_COLUMNS = {
		{Waves::T::P_FORWARD,   1},
		{Waves::T::P_BACKWARD,  0},
		{Waves::T::S1_FORWARD,  4},
		{Waves::T::S1_BACKWARD, 2},
		{Waves::T::S2_FORWARD,  5},
		{Waves::T::S2_BACKWARD, 3}
};

template<>
const std::map<Waves::T, int/* number of column in U1 */> GcmMatrices<VelocitySigmaVariables<3>, OrthotropicMaterial>::WAVE_COLUMNS = {
		{Waves::T::P_FORWARD,   5},
		{Waves::T::P_BACKWARD,  4},
		{Waves::T::S1_FORWARD,  1},
		{Waves::T::S1_BACKWARD, 0},
		{Waves::T::S2_FORWARD,  3},
		{Waves::T::S2_BACKWARD, 2}
};

template<int M>
real GcmMatrix<M>::getMaximalEigenvalue() const {
	real ans = 0.0;
	for (int i = 0; i < M; i++) {
		ans = fmax(ans, fabs(L(i, i)));
	}
	return ans;
};

template<typename TVariables, class TMaterial>
real GcmMatrices<TVariables, TMaterial>::getMaximalEigenvalue() const {
	real ans = 0.0;
	for (int i = 0; i < DIMENSIONALITY; i++) {
		ans = fmax(ans, A(i).getMaximalEigenvalue());
	}
	return ans;
};


template class GcmMatrices<VelocitySigmaVariables<1>, IsotropicMaterial>;
template class GcmMatrices<VelocitySigmaVariables<2>, IsotropicMaterial>;
template class GcmMatrices<VelocitySigmaVariables<3>, IsotropicMaterial>;
template class GcmMatrices<VelocitySigmaVariables<3>, OrthotropicMaterial>;