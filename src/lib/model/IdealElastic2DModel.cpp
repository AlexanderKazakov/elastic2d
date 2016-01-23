#include <lib/model/IdealElastic2DModel.hpp>

using namespace gcm;

const std::map<PhysicalQuantities::T, IdealElastic2DModel::GetSet> IdealElastic2DModel::QUANTITIES = {
		{PhysicalQuantities::T::Vx,       GetSet(IdealElastic2DModel::GetVx, IdealElastic2DModel::SetVx)},
		{PhysicalQuantities::T::Vy,       GetSet(IdealElastic2DModel::GetVy, IdealElastic2DModel::SetVy)},
		{PhysicalQuantities::T::Sxx,      GetSet(IdealElastic2DModel::GetSxx, IdealElastic2DModel::SetSxx)},
		{PhysicalQuantities::T::Sxy,      GetSet(IdealElastic2DModel::GetSxy, IdealElastic2DModel::SetSxy)},
		{PhysicalQuantities::T::Syy,      GetSet(IdealElastic2DModel::GetSyy, IdealElastic2DModel::SetSyy)},
		{PhysicalQuantities::T::PRESSURE, GetSet(IdealElastic2DModel::GetPressure, IdealElastic2DModel::SetPressure)}
};
