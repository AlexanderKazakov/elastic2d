#include "IdealElastic1DModel.hpp"

using namespace gcm;

const std::map<PhysicalQuantities::QUANTITY, IdealElastic1DModel::GetSet> IdealElastic1DModel::QUANTITIES = {
		{PhysicalQuantities::QUANTITY::Vx,       GetSet(IdealElastic1DModel::GetVx, IdealElastic1DModel::SetVx)},
		{PhysicalQuantities::QUANTITY::Sxx,      GetSet(IdealElastic1DModel::GetSxx, IdealElastic1DModel::SetSxx)},
		{PhysicalQuantities::QUANTITY::PRESSURE, GetSet(IdealElastic1DModel::GetPressure, IdealElastic1DModel::SetPressure)}
};
