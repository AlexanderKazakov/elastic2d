#ifndef LIBGCM_LINAL_VECTOR3_HPP
#define LIBGCM_LINAL_VECTOR3_HPP

#include <lib/linal/Linal.hpp>

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

		Vector<3, Vector3Container> crossProduct(const Vector<3, Vector3Container> &v1, const Vector<3, Vector3Container> &v2);

		typedef Vector<3, Vector3Container> Vector3;
	};
};

#endif // LIBGCM_LINAL_VECTOR3_HPP
