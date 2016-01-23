#include <lib/model/IdealElastic3DModel.hpp>

using namespace gcm;


const std::map<PhysicalQuantities::T, IdealElastic3DModel::GetSet> IdealElastic3DModel::QUANTITIES = {
		{PhysicalQuantities::T::Vx,       GetSet(IdealElastic3DModel::GetVx, IdealElastic3DModel::SetVx)},
		{PhysicalQuantities::T::Vy,       GetSet(IdealElastic3DModel::GetVy, IdealElastic3DModel::SetVy)},
		{PhysicalQuantities::T::Vz,       GetSet(IdealElastic3DModel::GetVz, IdealElastic3DModel::SetVz)},
		{PhysicalQuantities::T::Sxx,      GetSet(IdealElastic3DModel::GetSxx, IdealElastic3DModel::SetSxx)},
		{PhysicalQuantities::T::Sxy,      GetSet(IdealElastic3DModel::GetSxy, IdealElastic3DModel::SetSxy)},
		{PhysicalQuantities::T::Sxz,      GetSet(IdealElastic3DModel::GetSxz, IdealElastic3DModel::SetSxz)},
		{PhysicalQuantities::T::Syy,      GetSet(IdealElastic3DModel::GetSyy, IdealElastic3DModel::SetSyy)},
		{PhysicalQuantities::T::Syz,      GetSet(IdealElastic3DModel::GetSyz, IdealElastic3DModel::SetSyz)},
		{PhysicalQuantities::T::Szz,      GetSet(IdealElastic3DModel::GetSzz, IdealElastic3DModel::SetSzz)},
		{PhysicalQuantities::T::PRESSURE, GetSet(IdealElastic3DModel::GetPressure, IdealElastic3DModel::SetPressure)}
};
