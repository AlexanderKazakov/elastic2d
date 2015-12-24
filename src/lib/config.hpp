#ifndef LIBGCM_CONFIG_HPP
#define LIBGCM_CONFIG_HPP

#define LIBGCM_DOUBLE_PRECISION 1


#define CONFIG_ENABLE_ASSERTIONS 1

#define CONFIG_ENABLE_LOGGING 1
#define CONFIG_ENABLE_LOGGING_TRACE 1
#define CONFIG_ENABLE_LOGGING_ERROR 1
#define CONFIG_ENABLE_LOGGING_FATAL 1
#define CONFIG_ENABLE_LOGGING_WARN 1
#define CONFIG_ENABLE_LOGGING_INFO 1
#define CONFIG_ENABLE_LOGGING_DEBUG 1


const int N = 5; // number of variables in PDE
enum class NodeMap {Vx, Vy, Sxx, Sxy, Syy}; // order of physical quantities in the vector in Node

enum class InitialConditions {Zero, TestExplosion, Explosion, ExplosionAtTheLeft, PWaveXBackward, PWaveYBackward,
	PWaveX, SWaveX, PWaveY, SWaveY, SWaveXBackward, SxxOnly};

enum class BorderConditions {NonReflection, FreeBorder};


#endif //LIBGCM_CONFIG_HPP
