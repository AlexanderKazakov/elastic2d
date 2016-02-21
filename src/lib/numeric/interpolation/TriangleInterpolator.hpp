#ifndef LIBGCM_TRIANGLEINTERPOLATOR_HPP
#define LIBGCM_TRIANGLEINTERPOLATOR_HPP

#include <vector>
#include <lib/linal/linal.hpp>
#include <lib/util/Types.hpp>


namespace gcm {
	template<class TValue>
	class TriangleInterpolator {
		// TODO - some non-dimensionalization?
	public:
		typedef linal::Vector2 Coords;

		/**
		 * Linear interpolation in plain triangle
		 * @param c_i and v_i - points and values
		 * @param q point to interpolate
		 */
		static TValue interpolate(const Coords& c1, const TValue v1,
		                          const Coords& c2, const TValue v2,
		                          const Coords& c3, const TValue v3,
								  const Coords& q) {

			real det = linal::determinant(c1.x, c1.y, 1.0,
			                              c2.x, c2.y, 1.0,
										  c3.x, c3.y, 1.0);

			TValue det1 = linal::determinant(v1, c1.y, 1,
			                                 v2, c2.y, 1,
										     v3, c3.y, 1);

			TValue det2 = linal::determinant(c1.x, v1, 1,
			                                 c2.x, v2, 1,
										     c3.x, v3, 1);

			TValue det3 = linal::determinant(c1.x, c1.y, v1,
			                                 c2.x, c2.y, v2,
										     c3.x, c3.y, v3);

			// v(x, y) = ax + by + c
			auto a = det1 / det;
			auto b = det2 / det;
			auto c = det3 / det;

			return a * q.x + b * q.y + c;
		};
	};
}

#endif // LIBGCM_TRIANGLEINTERPOLATOR_HPP
