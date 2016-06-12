#include <lib/util/Enum.hpp>

using namespace gcm;
/* *INDENT-OFF* */


const std::map<PhysicalQuantities::T, std::string> PhysicalQuantities::NAME = {
		{PhysicalQuantities::T::VELOCITY,       "Velocity"},
		{PhysicalQuantities::T::FORCE,          "Force"},
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


const std::map<BorderConditions::T, std::string> BorderConditions::NAME = {
		{BorderConditions::T::FIXED_FORCE,    "fixed_force"},
		{BorderConditions::T::FIXED_VELOCITY, "fixed_velocity"}
};


const std::map<Materials::T, std::string> Materials::NAME = {
		{Materials::T::ISOTROPIC,    "isotropic"},
		{Materials::T::ORTHOTROPIC,  "orthotropic"}
};


const std::map<Grids::T, std::string> Grids::NAME = {
		{Grids::T::CUBIC,    "cubic"},
		{Grids::T::CGAL,     "cgal"}
};


const std::map<Models::T, std::string> Models::NAME = {
		{Models::T::ELASTIC1D,    "elastic1d"},
		{Models::T::ELASTIC2D,    "elastic2d"},
		{Models::T::ELASTIC3D,    "elastic3d"}
};


const std::map<Snapshotters::T, std::string> Snapshotters::NAME = {
		{Snapshotters::T::VTK,            "vtk"},
		{Snapshotters::T::BIN2DSEISM,     "bin2dseism"},
		{Snapshotters::T::DETECTOR,       "detector"},
		{Snapshotters::T::SLICESNAP,      "slice_snap"},
};

