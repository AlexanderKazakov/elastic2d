#include "IdealElastic3DModel.hpp"

using namespace gcm;


const std::map<PhysicalQuantities::QUANTITY, IdealElastic3DModel::GetSet> IdealElastic3DModel::QUANTITIES = {
		{PhysicalQuantities::QUANTITY::Vx,       GetSet(IdealElastic3DModel::GetVx, IdealElastic3DModel::SetVx)},
		{PhysicalQuantities::QUANTITY::Vy,       GetSet(IdealElastic3DModel::GetVy, IdealElastic3DModel::SetVy)},
		{PhysicalQuantities::QUANTITY::Vz,       GetSet(IdealElastic3DModel::GetVz, IdealElastic3DModel::SetVz)},
		{PhysicalQuantities::QUANTITY::Sxx,      GetSet(IdealElastic3DModel::GetSxx, IdealElastic3DModel::SetSxx)},
		{PhysicalQuantities::QUANTITY::Sxy,      GetSet(IdealElastic3DModel::GetSxy, IdealElastic3DModel::SetSxy)},
		{PhysicalQuantities::QUANTITY::Sxz,      GetSet(IdealElastic3DModel::GetSxz, IdealElastic3DModel::SetSxz)},
		{PhysicalQuantities::QUANTITY::Syy,      GetSet(IdealElastic3DModel::GetSyy, IdealElastic3DModel::SetSyy)},
		{PhysicalQuantities::QUANTITY::Syz,      GetSet(IdealElastic3DModel::GetSyz, IdealElastic3DModel::SetSyz)},
		{PhysicalQuantities::QUANTITY::Szz,      GetSet(IdealElastic3DModel::GetSzz, IdealElastic3DModel::SetSzz)},
		{PhysicalQuantities::QUANTITY::PRESSURE, GetSet(IdealElastic3DModel::GetPressure, IdealElastic3DModel::SetPressure)}
};
