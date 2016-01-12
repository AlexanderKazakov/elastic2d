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
#define CONFIG_ENABLE_LOGGING_DEBUG 0


enum class Border{X_LEFT = 0, X_RIGHT = 1, Y_LEFT = 2, Y_RIGHT = 3, Z_LEFT = 4, Z_RIGHT = 5};

enum class InitialConditions {Zero, TestExplosion, Explosion, PWaveXBackward, PWaveYBackward,
	PWaveX, SWaveX, PWaveY, SWaveY, SWaveXBackward, SxxOnly};
enum class BorderConditions {NonReflection, FreeBorder};


#endif //LIBGCM_CONFIG_HPP
