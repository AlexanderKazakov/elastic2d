#ifndef LIBGCM_ISOTROPICMATERIAL_HPP
#define LIBGCM_ISOTROPICMATERIAL_HPP

#include <lib/rheology/gcm_matrices/GcmMatrices.hpp>
#include <lib/rheology/variables/variables.hpp>

namespace gcm {
	class Task;

	class IsotropicMaterial {
	public:
		real rho = 0; // density

		real lambda = 0; // Lame parameters
		real mu = 0;

		real yieldStrength = 0; // plasticity parameter

		real continualDamageParameter = 0; // parameter in continual damage equation

		IsotropicMaterial(const IsotropicMaterial& other) = default;
		IsotropicMaterial(const real _rho = 0, const real _lambda = 0, const real _mu = 0,
		                  const real _yieldStrength = 0, const real _continualDamageParameter = 0);
		
		void initialize(const Task& task);

		/** Fill in gcm matrices */
		void constructGcmMatrices(GcmMatrices<VelocitySigmaVariables<1>, IsotropicMaterial>& m) const;
		void constructGcmMatrices(GcmMatrices<VelocitySigmaVariables<2>, IsotropicMaterial>& m) const;
		void constructGcmMatrices(GcmMatrices<VelocitySigmaVariables<3>, IsotropicMaterial>& m) const;

		/** For testing purposes */
		static IsotropicMaterial generateRandomMaterial();

	};
}

#endif // LIBGCM_ISOTROPICMATERIAL_HPP
