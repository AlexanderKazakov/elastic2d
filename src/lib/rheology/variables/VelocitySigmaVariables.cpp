#include <lib/rheology/variables/VelocitySigmaVariables.hpp>

using namespace gcm;

#define GETSETTER_VELOCITY_PAIR(q, d, i) {PhysicalQuantities::T::q, GetSetter<VelocitySigmaVariables<d>> \
		(VelocitySigmaVariables<d>::GetVelocity<i>, VelocitySigmaVariables<d>::SetVelocity<i>)}

#define GETSETTER_SIGMA_PAIR(q, d, i, j) {PhysicalQuantities::T::q, GetSetter<VelocitySigmaVariables<d>> \
		(VelocitySigmaVariables<d>::GetSigma<i, j>, VelocitySigmaVariables<d>::SetSigma<i, j>)}


template<> const std::map<PhysicalQuantities::T, GetSetter<VelocitySigmaVariables<1>>>
		VelocitySigmaVariables<1>::QUANTITIES = {

		GETSETTER_VELOCITY_PAIR(Vx, 1, 0),
        GETSETTER_SIGMA_PAIR(Sxx, 1, 0, 0),
		{PhysicalQuantities::T::PRESSURE, GetSetter<VelocitySigmaVariables<1>>
				(VelocitySigmaVariables<1>::GetPressure, VelocitySigmaVariables<1>::SetPressure)}
};

template<> const std::map<PhysicalQuantities::T, GetSetter<VelocitySigmaVariables<2>>>
		VelocitySigmaVariables<2>::QUANTITIES = {

		GETSETTER_VELOCITY_PAIR(Vx, 2, 0),
		GETSETTER_VELOCITY_PAIR(Vy, 2, 1),

		GETSETTER_SIGMA_PAIR(Sxx, 2, 0, 0),
		GETSETTER_SIGMA_PAIR(Sxy, 2, 0, 1),
		GETSETTER_SIGMA_PAIR(Syy, 2, 1, 1),
		{PhysicalQuantities::T::PRESSURE, GetSetter<VelocitySigmaVariables<2>>
				(VelocitySigmaVariables<2>::GetPressure, VelocitySigmaVariables<2>::SetPressure)}
};

template<> const std::map<PhysicalQuantities::T, GetSetter<VelocitySigmaVariables<3>>>
		VelocitySigmaVariables<3>::QUANTITIES = {

		GETSETTER_VELOCITY_PAIR(Vx, 3, 0),
		GETSETTER_VELOCITY_PAIR(Vy, 3, 1),
		GETSETTER_VELOCITY_PAIR(Vz, 3, 2),

		GETSETTER_SIGMA_PAIR(Sxx, 3, 0, 0),
		GETSETTER_SIGMA_PAIR(Sxy, 3, 0, 1),
		GETSETTER_SIGMA_PAIR(Sxz, 3, 0, 2),
		GETSETTER_SIGMA_PAIR(Syy, 3, 1, 1),
		GETSETTER_SIGMA_PAIR(Syz, 3, 1, 2),
		GETSETTER_SIGMA_PAIR(Szz, 3, 2, 2),
		{PhysicalQuantities::T::PRESSURE, GetSetter<VelocitySigmaVariables<3>>
				(VelocitySigmaVariables<3>::GetPressure, VelocitySigmaVariables<3>::SetPressure)}
};

#undef GETSETTER_SIGMA_PAIR
#undef GETSETTER_VELOCITY_PAIR

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