#include <lib/rheology/ode/ode.hpp>

using namespace gcm;

const std::map<PhysicalQuantities::T, GetSetter<ContinualDamageOde::Variables>> ContinualDamageOde::QUANTITIES = {
		{PhysicalQuantities::T::DAMAGE_MEASURE,
				GetSetter<ContinualDamageOde::Variables>(ContinualDamageOde::GetHi, ContinualDamageOde::SetHi)}
};
