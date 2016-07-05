#ifndef LIBGCM_MODEL_HPP
#define LIBGCM_MODEL_HPP

#include <lib/rheology/materials/materials.hpp>
#include <lib/rheology/variables/variables.hpp>
#include <lib/rheology/ode/ode.hpp>
#include <lib/rheology/correctors/correctors.hpp>
#include <lib/numeric/gcm/GcmMatrices.hpp>


namespace gcm {

/**
 * Map between wave types and corresponding columns in matrix U1 in GcmMatrix.
 */
typedef std::map<Waves::T, int>                      WavesEigenvectorsMap;


/**
 * Map between material types and corresponding WavesEigenvectorsMap for concrete Model.
 * Note that in concrete Model for concrete Material the order of eigenvalues
 * in GcmMatrices has to be the same for all spatial directions.
 */
typedef std::map<Materials::T, WavesEigenvectorsMap> MaterialsWavesMap;


}


#endif // LIBGCM_MODEL_HPP
