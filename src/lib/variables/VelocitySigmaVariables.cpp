#include <lib/variables/VelocitySigmaVariables.hpp>

using namespace gcm;

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

template class VelocitySigmaVariables<1>;
template class VelocitySigmaVariables<2>;
template class VelocitySigmaVariables<3>;