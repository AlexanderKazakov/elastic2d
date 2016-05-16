#ifndef LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
#define LIBGCM_GRIDCHARACTERISTICMETHOD_HPP

#include <lib/linal/linal.hpp>
#include <lib/numeric/interpolation/interpolation.hpp>
#include <lib/mesh/DefaultMesh.hpp>


namespace gcm {

/** @return values on next time layer by grid-characteristic method */
template<typename TMatrix>
inline auto
localGcmStep(const TMatrix& U1, const TMatrix& U, const TMatrix& values) 
		-> decltype(U1 * linal::diagonalMultiply(U, values)) {
	return /* new values = U1 * Riemann solvers */ U1 * linal::diagonalMultiply(
	       /* Riemann solvers = U * old values */ U,
	       /* old values are in columns of the matrix */ values);
}


/**
 * Grid-characteristic method
 */
template<typename TModel, typename TGrid, typename TMaterial>
class GridCharacteristicMethod;


}

#endif // LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
