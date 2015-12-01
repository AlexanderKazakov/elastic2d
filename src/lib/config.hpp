#ifndef ELASTIC2D_CONFIG_HPP
#define ELASTIC2D_CONFIG_HPP

#define ELASTIC2D_DOUBLE_PRECISION 1

#if ELASTIC2D_DOUBLE_PRECISION
typedef double real;
const real EQUALITY_TOLERANCE = 1e-9;
#else
typedef float real;
const real EQUALITY_TOLERANCE = 1e-4;
#endif

typedef unsigned int uint;

const uint N = 5; // number of variables in PDE
enum class NodeMap {Vx, Vy, Sxx, Sxy, Syy}; // order of physical quantities in the vector in Node

enum class InitialConditions {Zero, TestExplosion, Explosion, PWaveX, SWaveX, PWaveY, SWaveY,
	SWaveXBackward, SxxOnly};


#endif //ELASTIC2D_CONFIG_HPP
