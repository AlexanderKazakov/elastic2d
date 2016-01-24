#ifndef LIBGCM_LINAL_VECTOR2_HPP
#define LIBGCM_LINAL_VECTOR2_HPP

#include <lib/linal/Linal.hpp>

namespace gcm {
	namespace linal {
		/**
		 * Specialized value container for 2-dimensional vector.
		 */
		struct Vector2Container {
			static const int SIZE = 2; // size of storage in units of gcm::real
			union {
				real values[SIZE];
				struct {
					real x;
					real y;
				};
			};
		};


		/**
		 * Specialized 2-dimensional vector
		 */
		typedef Vector<2, Vector2Container> Vector2;
	};
};

#endif // LIBGCM_LINAL_VECTOR2_HPP
