#include <libgcm/rheology/variables/AcousticVariables.hpp>

using namespace gcm;


#define GETSETTER_V_PAIR(q, d, i) {PhysicalQuantities::T::q, \
                                   AcousticVariables<d>::GETSETTER \
                                       (AcousticVariables<d>::GetV<i>, \
                                        AcousticVariables<d>::SetV<i>)}

#define VELOCITY_GETSETTER(d) \
template<> const AcousticVariables<d>::Vector3Map \
		AcousticVariables<d>::VECTORS = { \
				{PhysicalQuantities::T::VELOCITY, AcousticVariables<d>::VECTOR3GETSETTER \
						(AcousticVariables<d>::GetVelocity, \
						 AcousticVariables<d>::SetVelocity)} }


template<> const AcousticVariables<1>::QuantitiesMap
		AcousticVariables<1>::QUANTITIES = {

		GETSETTER_V_PAIR(Vx, 1, 0),
		{PhysicalQuantities::T::PRESSURE, AcousticVariables<1>::GETSETTER
				(AcousticVariables<1>::GetPressure,
				 AcousticVariables<1>::SetPressure)}
};
VELOCITY_GETSETTER(1);

template<> const AcousticVariables<2>::QuantitiesMap
		AcousticVariables<2>::QUANTITIES = {

		GETSETTER_V_PAIR(Vx, 2, 0),
		GETSETTER_V_PAIR(Vy, 2, 1),

		{PhysicalQuantities::T::PRESSURE, AcousticVariables<2>::GETSETTER
				(AcousticVariables<2>::GetPressure,
				 AcousticVariables<2>::SetPressure)}
};
VELOCITY_GETSETTER(2);

template<> const AcousticVariables<3>::QuantitiesMap
		AcousticVariables<3>::QUANTITIES = {

		GETSETTER_V_PAIR(Vx, 3, 0),
		GETSETTER_V_PAIR(Vy, 3, 1),
		GETSETTER_V_PAIR(Vz, 3, 2),

		{PhysicalQuantities::T::PRESSURE, AcousticVariables<3>::GETSETTER
				(AcousticVariables<3>::GetPressure,
				 AcousticVariables<3>::SetPressure)}
};
VELOCITY_GETSETTER(3);

#undef GETSETTER_SIGMA_PAIR
#undef GETSETTER_V_PAIR
#undef VELOCITY_GETSETTER


template class AcousticVariables<1>;
template class AcousticVariables<2>;
template class AcousticVariables<3>;
