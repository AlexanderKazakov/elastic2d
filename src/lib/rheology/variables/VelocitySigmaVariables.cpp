#include <lib/rheology/variables/VelocitySigmaVariables.hpp>

using namespace gcm;

template<> const std::map<PhysicalQuantities::T, GetSetter<VelocitySigmaVariables<1>>> VelocitySigmaVariables<1>::QUANTITIES = {
		{PhysicalQuantities::T::Vx,       GetSetter<VelocitySigmaVariables<1>>(VelocitySigmaVariables<1>::GetVelocity<0>, VelocitySigmaVariables<1>::SetVelocity<0>)},
		{PhysicalQuantities::T::Sxx,      GetSetter<VelocitySigmaVariables<1>>(VelocitySigmaVariables<1>::GetSigma<0, 0>, VelocitySigmaVariables<1>::SetSigma<0, 0>)},
		{PhysicalQuantities::T::PRESSURE, GetSetter<VelocitySigmaVariables<1>>(VelocitySigmaVariables<1>::GetPressure, VelocitySigmaVariables<1>::SetPressure)}
};

template<> const std::map<PhysicalQuantities::T, GetSetter<VelocitySigmaVariables<2>>> VelocitySigmaVariables<2>::QUANTITIES = {
		{PhysicalQuantities::T::Vx,       GetSetter<VelocitySigmaVariables<2>>(VelocitySigmaVariables<2>::GetVelocity<0>, VelocitySigmaVariables<2>::SetVelocity<0>)},
		{PhysicalQuantities::T::Vy,       GetSetter<VelocitySigmaVariables<2>>(VelocitySigmaVariables<2>::GetVelocity<1>, VelocitySigmaVariables<2>::SetVelocity<1>)},
		{PhysicalQuantities::T::Sxx,      GetSetter<VelocitySigmaVariables<2>>(VelocitySigmaVariables<2>::GetSigma<0, 0>, VelocitySigmaVariables<2>::SetSigma<0, 0>)},
		{PhysicalQuantities::T::Sxy,      GetSetter<VelocitySigmaVariables<2>>(VelocitySigmaVariables<2>::GetSigma<0, 1>, VelocitySigmaVariables<2>::SetSigma<0, 1>)},
		{PhysicalQuantities::T::Syy,      GetSetter<VelocitySigmaVariables<2>>(VelocitySigmaVariables<2>::GetSigma<1, 1>, VelocitySigmaVariables<2>::SetSigma<1, 1>)},
		{PhysicalQuantities::T::PRESSURE, GetSetter<VelocitySigmaVariables<2>>(VelocitySigmaVariables<2>::GetPressure, VelocitySigmaVariables<2>::SetPressure)}
};

template<> const std::map<PhysicalQuantities::T, GetSetter<VelocitySigmaVariables<3>>> VelocitySigmaVariables<3>::QUANTITIES = {
		{PhysicalQuantities::T::Vx,       GetSetter<VelocitySigmaVariables<3>>(VelocitySigmaVariables<3>::GetVelocity<0>, VelocitySigmaVariables<3>::SetVelocity<0>)},
		{PhysicalQuantities::T::Vy,       GetSetter<VelocitySigmaVariables<3>>(VelocitySigmaVariables<3>::GetVelocity<1>, VelocitySigmaVariables<3>::SetVelocity<1>)},
		{PhysicalQuantities::T::Vz,       GetSetter<VelocitySigmaVariables<3>>(VelocitySigmaVariables<3>::GetVelocity<2>, VelocitySigmaVariables<3>::SetVelocity<2>)},
		{PhysicalQuantities::T::Sxx,      GetSetter<VelocitySigmaVariables<3>>(VelocitySigmaVariables<3>::GetSigma<0, 0>, VelocitySigmaVariables<3>::SetSigma<0, 0>)},
		{PhysicalQuantities::T::Sxy,      GetSetter<VelocitySigmaVariables<3>>(VelocitySigmaVariables<3>::GetSigma<0, 1>, VelocitySigmaVariables<3>::SetSigma<0, 1>)},
		{PhysicalQuantities::T::Sxz,      GetSetter<VelocitySigmaVariables<3>>(VelocitySigmaVariables<3>::GetSigma<0, 2>, VelocitySigmaVariables<3>::SetSigma<0, 2>)},
		{PhysicalQuantities::T::Syy,      GetSetter<VelocitySigmaVariables<3>>(VelocitySigmaVariables<3>::GetSigma<1, 1>, VelocitySigmaVariables<3>::SetSigma<1, 1>)},
		{PhysicalQuantities::T::Syz,      GetSetter<VelocitySigmaVariables<3>>(VelocitySigmaVariables<3>::GetSigma<1, 2>, VelocitySigmaVariables<3>::SetSigma<1, 2>)},
		{PhysicalQuantities::T::Szz,      GetSetter<VelocitySigmaVariables<3>>(VelocitySigmaVariables<3>::GetSigma<2, 2>, VelocitySigmaVariables<3>::SetSigma<2, 2>)},
		{PhysicalQuantities::T::PRESSURE, GetSetter<VelocitySigmaVariables<3>>(VelocitySigmaVariables<3>::GetPressure, VelocitySigmaVariables<3>::SetPressure)}
};


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