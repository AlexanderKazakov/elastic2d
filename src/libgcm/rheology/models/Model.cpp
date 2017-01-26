#include <libgcm/rheology/models/models.hpp>

using namespace gcm;


template<>
const MaterialsWavesMap ElasticModel<1>::MATERIALS_WAVES_MAP = {
		{(Materials::T) IsotropicMaterial::Type, {
				{Waves::T::P_FORWARD,  0},
				{Waves::T::P_BACKWARD, 1}
		 }}
};


template<>
const MaterialsWavesMap ElasticModel<2>::MATERIALS_WAVES_MAP = {
		{(Materials::T) IsotropicMaterial::Type,	{
				{Waves::T::P_FORWARD,   0},
				{Waves::T::P_BACKWARD,  1},
				{Waves::T::S1_FORWARD,  2},
				{Waves::T::S1_BACKWARD, 3}
		 }},
		{(Materials::T) OrthotropicMaterial::Type, {
				 {Waves::T::P_FORWARD,   3},
				 {Waves::T::P_BACKWARD,  2},
				 {Waves::T::S1_FORWARD,  1},
				 {Waves::T::S1_BACKWARD, 0},
		 }}
};


template<>
const MaterialsWavesMap ElasticModel<3>::MATERIALS_WAVES_MAP = {
		{(Materials::T) IsotropicMaterial::Type, {
				{Waves::T::P_FORWARD,   0},
				{Waves::T::P_BACKWARD,  1},
				{Waves::T::S1_FORWARD,  2},
				{Waves::T::S1_BACKWARD, 3},
				{Waves::T::S2_FORWARD,  4},
				{Waves::T::S2_BACKWARD, 5}
		 }},
		{(Materials::T) OrthotropicMaterial::Type, {
				 {Waves::T::P_FORWARD,   5},
				 {Waves::T::P_BACKWARD,  4},
				 {Waves::T::S1_FORWARD,  1},
				 {Waves::T::S1_BACKWARD, 0},
				 {Waves::T::S2_FORWARD,  3},
				 {Waves::T::S2_BACKWARD, 2}
		 }}
};



template<>
const MaterialsWavesMap AcousticModel<3>::MATERIALS_WAVES_MAP = {
		{(Materials::T) IsotropicMaterial::Type, {
				{Waves::T::P_FORWARD,   0},
				{Waves::T::P_BACKWARD,  1},
		 }}
};

template<>
const MaterialsWavesMap AcousticModel<2>::MATERIALS_WAVES_MAP = 
		AcousticModel<3>::MATERIALS_WAVES_MAP;
template<>
const MaterialsWavesMap AcousticModel<1>::MATERIALS_WAVES_MAP = 
		AcousticModel<3>::MATERIALS_WAVES_MAP;


template<> const AcousticModel<1>::WaveIndices AcousticModel<1>:: POSITIVE_WAVES = {0};
template<> const AcousticModel<1>::WaveIndices AcousticModel<1>::NEGATIVE_WAVES = {1};
template<> const AcousticModel<2>::WaveIndices AcousticModel<2>:: POSITIVE_WAVES = {0};
template<> const AcousticModel<2>::WaveIndices AcousticModel<2>::NEGATIVE_WAVES = {1};
template<> const AcousticModel<3>::WaveIndices AcousticModel<3>:: POSITIVE_WAVES = {0};
template<> const AcousticModel<3>::WaveIndices AcousticModel<3>::NEGATIVE_WAVES = {1};

template<> const ElasticModel<1>::WaveIndices ElasticModel<1>:: POSITIVE_WAVES = {0};
template<> const ElasticModel<1>::WaveIndices ElasticModel<1>::NEGATIVE_WAVES = {1};
template<> const ElasticModel<2>::WaveIndices ElasticModel<2>:: POSITIVE_WAVES = {0, 2};
template<> const ElasticModel<2>::WaveIndices ElasticModel<2>::NEGATIVE_WAVES = {1, 3};
template<> const ElasticModel<3>::WaveIndices ElasticModel<3>:: POSITIVE_WAVES = {0, 2, 4};
template<> const ElasticModel<3>::WaveIndices ElasticModel<3>::NEGATIVE_WAVES = {1, 3, 5};

