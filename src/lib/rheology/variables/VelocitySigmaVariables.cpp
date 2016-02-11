#include <lib/rheology/variables/VelocitySigmaVariables.hpp>

using namespace gcm;

#define GETSETTER_V_PAIR(q, d, i) {PhysicalQuantities::T::q, GetSetter<VelocitySigmaVariables<d>> \
		(VelocitySigmaVariables<d>::GetV<i>, VelocitySigmaVariables<d>::SetV<i>)}

#define GETSETTER_SIGMA_PAIR(q, d, i, j) {PhysicalQuantities::T::q, GetSetter<VelocitySigmaVariables<d>> \
		(VelocitySigmaVariables<d>::GetSigma<i, j>, VelocitySigmaVariables<d>::SetSigma<i, j>)}

#define VELOCITY_GETSETTER(d) \
template<> const std::map<PhysicalQuantities::T, Vector3GetSetter<VelocitySigmaVariables<d>>> \
		VelocitySigmaVariables<d>::VECTORS = { \
				{PhysicalQuantities::T::VELOCITY, Vector3GetSetter<VelocitySigmaVariables<d>> \
						(VelocitySigmaVariables<d>::GetVelocity, VelocitySigmaVariables<d>::SetVelocity)} };


template<> const std::map<PhysicalQuantities::T, GetSetter<VelocitySigmaVariables<1>>>
		VelocitySigmaVariables<1>::QUANTITIES = {

		GETSETTER_V_PAIR(Vx, 1, 0),
        GETSETTER_SIGMA_PAIR(Sxx, 1, 0, 0),
		{PhysicalQuantities::T::PRESSURE, GetSetter<VelocitySigmaVariables<1>>
				(VelocitySigmaVariables<1>::GetPressure, VelocitySigmaVariables<1>::SetPressure)}
};
VELOCITY_GETSETTER(1)

template<> const std::map<PhysicalQuantities::T, GetSetter<VelocitySigmaVariables<2>>>
		VelocitySigmaVariables<2>::QUANTITIES = {

		GETSETTER_V_PAIR(Vx, 2, 0),
		GETSETTER_V_PAIR(Vy, 2, 1),

		GETSETTER_SIGMA_PAIR(Sxx, 2, 0, 0),
		GETSETTER_SIGMA_PAIR(Sxy, 2, 0, 1),
		GETSETTER_SIGMA_PAIR(Syy, 2, 1, 1),
		{PhysicalQuantities::T::PRESSURE, GetSetter<VelocitySigmaVariables<2>>
				(VelocitySigmaVariables<2>::GetPressure, VelocitySigmaVariables<2>::SetPressure)}
};
VELOCITY_GETSETTER(2)

template<> const std::map<PhysicalQuantities::T, GetSetter<VelocitySigmaVariables<3>>>
		VelocitySigmaVariables<3>::QUANTITIES = {

		GETSETTER_V_PAIR(Vx, 3, 0),
		GETSETTER_V_PAIR(Vy, 3, 1),
		GETSETTER_V_PAIR(Vz, 3, 2),

		GETSETTER_SIGMA_PAIR(Sxx, 3, 0, 0),
		GETSETTER_SIGMA_PAIR(Sxy, 3, 0, 1),
		GETSETTER_SIGMA_PAIR(Sxz, 3, 0, 2),
		GETSETTER_SIGMA_PAIR(Syy, 3, 1, 1),
		GETSETTER_SIGMA_PAIR(Syz, 3, 1, 2),
		GETSETTER_SIGMA_PAIR(Szz, 3, 2, 2),
		{PhysicalQuantities::T::PRESSURE, GetSetter<VelocitySigmaVariables<3>>
				(VelocitySigmaVariables<3>::GetPressure, VelocitySigmaVariables<3>::SetPressure)}
};
VELOCITY_GETSETTER(3)

#undef GETSETTER_SIGMA_PAIR
#undef GETSETTER_V_PAIR
#undef VELOCITY_GETSETTER


template<int Dimensionality>
real VelocitySigmaVariables<Dimensionality>::getPressure() const {
	real trace = 0.0;
	for (int i = 0; i < Dimensionality; i++) {
		trace += this->sigma(i, i);
	}
	return - trace / Dimensionality;
}

template<int Dimensionality>
void VelocitySigmaVariables<Dimensionality>::setPressure(const real &pressure) {
	linal::clear(*this);
	for (int i = 0; i < Dimensionality; i++) {
		sigma(i, i) = - pressure;
	}
}

template<int Dimensionality>
real VelocitySigmaVariables<Dimensionality>::getJ2() const {
	real J22 = 0;
	real pressure = getPressure();
	for (int i = 0; i < Dimensionality; i++) {
		for (int j = 0; j < Dimensionality; j++) {
			J22 += (sigma(i, j) + (i == j) * pressure) * (sigma(i, j) + (i == j) * pressure) / 2;
		}
	}
	return sqrt(J22);
}

template class VelocitySigmaVariables<1>;
template class VelocitySigmaVariables<2>;
template class VelocitySigmaVariables<3>;