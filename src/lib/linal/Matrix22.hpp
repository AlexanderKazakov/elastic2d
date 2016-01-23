#ifndef LIBGCM_LINAL_MATRIX22_HPP
#define LIBGCM_LINAL_MATRIX22_HPP

#include <lib/linal/Linal.hpp>

namespace gcm {
	namespace linal {
		/**
		 * Specialized value container for 2x2 matrix.
		 */
		struct Matrix22Container {
			static const int SIZE = 2 * 2; // size of storage in units of gcm::real
			union {
				real values[SIZE];
				struct {
					real xx;
					real xy;
					real yx;
					real yy;
				};
			};
		};


		/**
		 * Specialized 2x2 matrix implementation.
		 */
		typedef Matrix<2, 2, Matrix22Container> Matrix22;

		/**
		 * @return matrix determinant
		 */
		real determinant(const Matrix22 &m);
	};
};

#endif // LIBGCM_LINAL_MATRIX22_HPP
