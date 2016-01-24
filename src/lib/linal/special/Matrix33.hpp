#ifndef LIBGCM_LINAL_MATRIX33_HPP
#define LIBGCM_LINAL_MATRIX33_HPP

#include <lib/linal/Linal.hpp>

namespace gcm {
	namespace linal {
		/**
		 * Specialized value container for 3x3 matrix.
		 */
		class Matrix33Container {
		public:
			static const int SIZE = 3 * 3; // size of storage in units of gcm::real
			union {
				real values[SIZE];
				struct {
					union {
						real xx;
						real a11;
					};
					union {
						real xy;
						real a12;
					};
					union {
						real xz;
						real a13;
					};
					union {
						real yx;
						real a21;
					};
					union {
						real yy;
						real a22;
					};
					union {
						real yz;
						real a23;
					};
					union {
						real zx;
						real a31;
					};
					union {
						real zy;
						real a32;
					};
					union {
						real zz;
						real a33;
					};
				};
			};
		};

		/**
		 * Specialized 3x3 matrix implementation.
		 */
		typedef Matrix<3, 3, Matrix33Container> Matrix33;

		/**
		 * @return matrix determinant
		 */
		real determinant(const Matrix33 &m);
	};
};

#endif // LIBGCM_LINAL_MATRIX33_HPP
