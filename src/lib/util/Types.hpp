#ifndef LIBGCM_TYPES_HPP
#define LIBGCM_TYPES_HPP

#include <lib/config.hpp>

namespace gcm {

#if LIBGCM_DOUBLE_PRECISION
typedef double real;
constexpr real EQUALITY_TOLERANCE = 1e-9;
#else
typedef float  real;
constexpr real EQUALITY_TOLERANCE = 1e-3;
#endif

}

#endif // LIBGCM_TYPES_HPP
