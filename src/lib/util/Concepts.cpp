#include <lib/util/Concepts.hpp>

using namespace gcm;


const std::map<PhysicalQuantities::T, std::string> PhysicalQuantities::NAME = {
		{PhysicalQuantities::T::Vx,             "Vx"},
		{PhysicalQuantities::T::Vy,             "Vy"},
		{PhysicalQuantities::T::Vz,             "Vz"},

		{PhysicalQuantities::T::Sxx,            "Sxx"},
		{PhysicalQuantities::T::Sxy,            "Sxy"},
		{PhysicalQuantities::T::Sxz,            "Sxz"},
		{PhysicalQuantities::T::Syy,            "Syy"},
		{PhysicalQuantities::T::Syz,            "Syz"},
		{PhysicalQuantities::T::Szz,            "Szz"},

		{PhysicalQuantities::T::RHO,            "rho"},
		{PhysicalQuantities::T::PRESSURE,       "pressure"},

		{PhysicalQuantities::T::DAMAGE_MEASURE, "damage_measure"},
};


const std::map<Waves::T, std::string> Waves::NAME = {
		{Waves::T::P_FORWARD,   "P_wave_forward"},
		{Waves::T::P_BACKWARD,  "P_wave_backward"},

		{Waves::T::S1_FORWARD,  "S1_wave_forward"},
		{Waves::T::S1_BACKWARD, "S1_wave_backward"},

		{Waves::T::S2_FORWARD,  "S2_wave_forward"},
		{Waves::T::S2_BACKWARD, "S2_wave_backward"},
};


const std::map<BorderCondition::T, std::string> BorderCondition::NAME = {
		{BorderCondition::T::NON_REFLECTION, "non-reflection"},
		{BorderCondition::T::FREE_BORDER,    "free-border"}
};

