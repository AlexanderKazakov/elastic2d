#ifndef LIBGCM_CONFIG_HPP
#define LIBGCM_CONFIG_HPP

#define LIBGCM_DOUBLE_PRECISION 1


#define CONFIG_ENABLE_ASSERTIONS 1

#define CONFIG_ENABLE_LOGGING 0

#define CONFIG_ENABLE_LOGGING_TRACE 1
#define CONFIG_ENABLE_LOGGING_ERROR 1
#define CONFIG_ENABLE_LOGGING_FATAL 1
#define CONFIG_ENABLE_LOGGING_WARN 1
#define CONFIG_ENABLE_LOGGING_INFO 1
#define CONFIG_ENABLE_LOGGING_DEBUG 1



enum class InitialConditions {Zero, TestExplosion, Explosion, ExplosionAtTheLeft, PWaveXBackward, PWaveYBackward,
	PWaveX, SWaveX, PWaveY, SWaveY, SWaveXBackward, SxxOnly};
enum class BorderConditions {NonReflection, FreeBorder};


#endif //LIBGCM_CONFIG_HPP
