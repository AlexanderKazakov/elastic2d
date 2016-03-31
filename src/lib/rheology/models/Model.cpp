#include <lib/rheology/models/Model.hpp>

using namespace gcm;
/* *INDENT-OFF* (disable code formatter) */

const MaterialsWavesMap Elastic1DModel::MATERIALS_WAVES_MAP = {
		{IsotropicMaterial::ID, {
				{Waves::T::P_FORWARD,  0},
				{Waves::T::P_BACKWARD, 1}
		 }}
};

const MaterialsWavesMap Elastic2DModel::MATERIALS_WAVES_MAP = {
		{IsotropicMaterial::ID,	{
				{Waves::T::P_FORWARD,   1},
				{Waves::T::P_BACKWARD,  0},
				{Waves::T::S1_FORWARD,  3},
				{Waves::T::S1_BACKWARD, 2}
		 }}
};

const MaterialsWavesMap Elastic3DModel::MATERIALS_WAVES_MAP = {
		{IsotropicMaterial::ID, {
				{Waves::T::P_FORWARD,   1},
				{Waves::T::P_BACKWARD,  0},
				{Waves::T::S1_FORWARD,  4},
				{Waves::T::S1_BACKWARD, 2},
				{Waves::T::S2_FORWARD,  5},
				{Waves::T::S2_BACKWARD, 3}
		 }},
		{OrthotropicMaterial::ID, {
				 {Waves::T::P_FORWARD,   5},
				 {Waves::T::P_BACKWARD,  4},
				 {Waves::T::S1_FORWARD,  1},
				 {Waves::T::S1_BACKWARD, 0},
				 {Waves::T::S2_FORWARD,  3},
				 {Waves::T::S2_BACKWARD, 2}
		 }}
};

const MaterialsWavesMap SuperDuperModel::MATERIALS_WAVES_MAP =
		Elastic3DModel::MATERIALS_WAVES_MAP;

