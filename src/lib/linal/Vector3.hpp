#ifndef LIBGCM_LINAL_VECTOR3_HPP
#define LIBGCM_LINAL_VECTOR3_HPP

#include "lib/linal/Linal.hpp"

namespace gcm {
	namespace linal {
		class Vector3Container {
		public:
			static const int SIZE = 3; // size of storage in units of gcm::real
			union {
				real values[SIZE];
				struct {
					union {
						real a1;
						real x;
					};
					union {
						real a2;
						real y;
					};
					union {
						real a3;
						real z;
					};
				};
			};
		};

		typedef Vector<3, Vector3Container> Vector3;

		Vector3 crossProduct(const Vector3 &v1, const Vector3 &v2);
	};
};

#endif // LIBGCM_LINAL_VECTOR3_HPP
