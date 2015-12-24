#ifndef LIBGCM_TYPES_HPP
#define LIBGCM_TYPES_HPP

#include "lib/config.hpp"

namespace gcm
{

#if LIBGCM_DOUBLE_PRECISION
	typedef double real;
	const real EQUALITY_TOLERANCE = 1e-9;
#else
	typedef float real;
	const real EQUALITY_TOLERANCE = 1e-4;
#endif

}

#endif // LIBGCM_TYPES_HPP
