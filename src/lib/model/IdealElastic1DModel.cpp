#include <lib/model/IdealElastic1DModel.hpp>

using namespace gcm;

const std::map<PhysicalQuantities::T, IdealElastic1DModel::GetSet> IdealElastic1DModel::QUANTITIES = {
		{PhysicalQuantities::T::Vx,       GetSet(IdealElastic1DModel::GetVx, IdealElastic1DModel::SetVx)},
		{PhysicalQuantities::T::Sxx,      GetSet(IdealElastic1DModel::GetSxx, IdealElastic1DModel::SetSxx)},
		{PhysicalQuantities::T::PRESSURE, GetSet(IdealElastic1DModel::GetPressure, IdealElastic1DModel::SetPressure)}
};
