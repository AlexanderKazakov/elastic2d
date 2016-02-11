#ifndef LIBGCM_COORDINATES_HPP
#define LIBGCM_COORDINATES_HPP

#include <lib/linal/linal.hpp>

namespace gcm {

	template<int Dimensionality>
	using CartesianCoords = linal::Vector<Dimensionality>;


	template<int Dimensionality>
	using DefaultCoords = CartesianCoords<Dimensionality>;
}

#endif // LIBGCM_COORDINATES_HPP
