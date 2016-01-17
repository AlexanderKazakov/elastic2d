#include "lib/util/Concepts.hpp"

using namespace gcm;


const std::map<PhysicalQuantities::QUANTITY, std::string> PhysicalQuantities::NAME = {
		{PhysicalQuantities::QUANTITY::Vx,       "Vx"},
		{PhysicalQuantities::QUANTITY::Vy,       "Vy"},
		{PhysicalQuantities::QUANTITY::Vz,       "Vz"},

		{PhysicalQuantities::QUANTITY::Sxx,      "Sxx"},
		{PhysicalQuantities::QUANTITY::Sxy,      "Sxy"},
		{PhysicalQuantities::QUANTITY::Sxz,      "Sxz"},
		{PhysicalQuantities::QUANTITY::Syy,      "Syy"},
		{PhysicalQuantities::QUANTITY::Syz,      "Syz"},
		{PhysicalQuantities::QUANTITY::Szz,      "Szz"},

		{PhysicalQuantities::QUANTITY::RHO,      "rho"},
		{PhysicalQuantities::QUANTITY::PRESSURE, "pressure"},
};

// TODO gtest
//static_assert(PhysicalQuantities::NAME.size() == static_cast<int>(PhysicalQuantities::QUANTITY::SIZE_OF_ENUM), "");


const std::map<Waves::WAVE, std::string> Waves::NAME = {
		{Waves::WAVE::P_FORWARD,   "P_wave_forward"},
		{Waves::WAVE::P_BACKWARD,  "P_wave_backward"},

		{Waves::WAVE::S1_FORWARD,  "S1_wave_forward"},
		{Waves::WAVE::S1_BACKWARD, "S1_wave_backward"},

		{Waves::WAVE::S2_FORWARD,  "S2_wave_forward"},
		{Waves::WAVE::S2_BACKWARD, "S2_wave_backward"},
};

//static_assert(Waves::NAME.size() == static_cast<int>(Waves::WAVE::SIZE_OF_ENUM), "");


const std::map<BorderCondition::CONDITION, std::string> BorderCondition::NAME = {
		{BorderCondition::CONDITION::NON_REFLECTION, "non-reflection"},
		{BorderCondition::CONDITION::FREE_BORDER,    "free-border"}
};

//static_assert(BorderCondition::NAME.size() == static_cast<int>(BorderCondition::CONDITION::SIZE_OF_ENUM), "");