#include "IdealElastic2DModel.hpp"

using namespace gcm;

const std::map<PhysicalQuantities::QUANTITY, IdealElastic2DModel::GetSet> IdealElastic2DModel::QUANTITIES = {
		{PhysicalQuantities::QUANTITY::Vx,       GetSet(IdealElastic2DModel::GetVx, IdealElastic2DModel::SetVx)},
		{PhysicalQuantities::QUANTITY::Vy,       GetSet(IdealElastic2DModel::GetVy, IdealElastic2DModel::SetVy)},
		{PhysicalQuantities::QUANTITY::Sxx,      GetSet(IdealElastic2DModel::GetSxx, IdealElastic2DModel::SetSxx)},
		{PhysicalQuantities::QUANTITY::Sxy,      GetSet(IdealElastic2DModel::GetSxy, IdealElastic2DModel::SetSxy)},
		{PhysicalQuantities::QUANTITY::Syy,      GetSet(IdealElastic2DModel::GetSyy, IdealElastic2DModel::SetSyy)},
		{PhysicalQuantities::QUANTITY::PRESSURE, GetSet(IdealElastic2DModel::GetPressure, IdealElastic2DModel::SetPressure)}
};
