#include "lib/linal/Vector3.hpp"

using namespace gcm;
using namespace gcm::linal;

Vector3 crossProduct(const Vector3 &v1, const Vector3 &v2) {
	return Vector3({v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x});
};