#ifndef LIBGCM_VELOCITYSIGMAVARIABLES_HPP
#define LIBGCM_VELOCITYSIGMAVARIABLES_HPP

#include <lib/linal/Linal.hpp>

namespace gcm {

	/**
	 * The most common gcm variables - velocity and symmetric tension tensor components
	 * @tparam Dimensionality space dimensionality
	 */
	template<int Dimensionality>
	struct VelocitySigmaVariables {
		static const int DIMENSIONALITY = Dimensionality;
		static const int SIZE = DIMENSIONALITY + ( DIMENSIONALITY * (DIMENSIONALITY + 1) ) / 2;
		union {
			real values[SIZE];
			struct {
				/** Velocity */
				real V[DIMENSIONALITY];
				/** Symmetric tension tensor components */
				real S[( DIMENSIONALITY * (DIMENSIONALITY + 1) ) / 2];
			};
		};

		/** Access to sigma as symmetric matrix */
		real sigma(const int i, const int j) const {
			return S[linal::SymmetricMatrix<Dimensionality>::getIndex(i, j)];
		};
		real& sigma(const int i, const int j) {
			return S[linal::SymmetricMatrix<Dimensionality>::getIndex(i, j)];
		};

		real getPressure() const;
		void setPressure(const real& pressure);
	};

}


#endif // LIBGCM_VELOCITYSIGMAVARIABLES_HPP
