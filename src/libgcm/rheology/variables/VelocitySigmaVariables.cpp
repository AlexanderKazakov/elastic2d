#include <libgcm/rheology/variables/VelocitySigmaVariables.hpp>

using namespace gcm;
/* *INDENT-OFF* */

#define GETSETTER_V_PAIR(q, d, i) {PhysicalQuantities::T::q, \
                                   VelocitySigmaVariables<d>::GETSETTER \
							           (VelocitySigmaVariables<d>::GetV<i>, \
                                        VelocitySigmaVariables<d>::SetV<i>)}

#define GETSETTER_SIGMA_PAIR(q, d, i, j) {PhysicalQuantities::T::q, \
                                          VelocitySigmaVariables<d>::GETSETTER \
                                              (VelocitySigmaVariables<d>::GetSigma<i, j>, \
                                               VelocitySigmaVariables<d>::SetSigma<i, j>)}

#define VELOCITY_GETSETTER(d) \
template<> const VelocitySigmaVariables<d>::Vector3Map \
		VelocitySigmaVariables<d>::VECTORS = { \
				{PhysicalQuantities::T::VELOCITY, VelocitySigmaVariables<d>::VECTOR3GETSETTER \
						(VelocitySigmaVariables<d>::GetVelocity, \
						 VelocitySigmaVariables<d>::SetVelocity)} }


template<> const VelocitySigmaVariables<1>::QuantitiesMap
		VelocitySigmaVariables<1>::QUANTITIES = {

		GETSETTER_V_PAIR(Vx, 1, 0),
        GETSETTER_SIGMA_PAIR(Sxx, 1, 0, 0),
		{PhysicalQuantities::T::PRESSURE, VelocitySigmaVariables<1>::GETSETTER
				(VelocitySigmaVariables<1>::GetPressure,
				 VelocitySigmaVariables<1>::SetPressure)}
};
VELOCITY_GETSETTER(1);

template<> const VelocitySigmaVariables<2>::QuantitiesMap
		VelocitySigmaVariables<2>::QUANTITIES = {

		GETSETTER_V_PAIR(Vx, 2, 0),
		GETSETTER_V_PAIR(Vy, 2, 1),

		GETSETTER_SIGMA_PAIR(Sxx, 2, 0, 0),
		GETSETTER_SIGMA_PAIR(Sxy, 2, 0, 1),
		GETSETTER_SIGMA_PAIR(Syy, 2, 1, 1),
		{PhysicalQuantities::T::PRESSURE, VelocitySigmaVariables<2>::GETSETTER
				(VelocitySigmaVariables<2>::GetPressure,
				 VelocitySigmaVariables<2>::SetPressure)}
};
VELOCITY_GETSETTER(2);

template<> const VelocitySigmaVariables<3>::QuantitiesMap
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
		{PhysicalQuantities::T::PRESSURE, VelocitySigmaVariables<3>::GETSETTER
				(VelocitySigmaVariables<3>::GetPressure,
				 VelocitySigmaVariables<3>::SetPressure)}
};
VELOCITY_GETSETTER(3);

#undef GETSETTER_SIGMA_PAIR
#undef GETSETTER_V_PAIR
#undef VELOCITY_GETSETTER


template class VelocitySigmaVariables<1>;
template class VelocitySigmaVariables<2>;
template class VelocitySigmaVariables<3>;
