#ifndef LIBGCM_ORTHOTROPICMATERIAL_HPP
#define LIBGCM_ORTHOTROPICMATERIAL_HPP

#include <initializer_list>

#include <lib/util/Types.hpp>
#include <lib/rheology/variables/variables.hpp>

namespace gcm {
	struct Statement;
	class IsotropicMaterial;

	class OrthotropicMaterial {
	public:
		static const Materials::T ID;

		real rho = 0; // density
		union {
			real c[9]; // elastic coefficients
			struct {
				real c11, c12, c13,
				          c22, c23,
				               c33,
				                    c44,
				                         c55,
				                              c66;
			};
		};

		real yieldStrength = 0; // plasticity parameters

		real continualDamageParameter = 0; // parameter in continual damage equation

		OrthotropicMaterial(const OrthotropicMaterial& other) = default;
		OrthotropicMaterial(const IsotropicMaterial& isotropic);
		OrthotropicMaterial(const real rho_ = 0, std::initializer_list<real> = {0,0,0,0,0,0,0,0,0},
		                    const real yieldStrength_ = 0, const real continualDamageParameter_ = 0);

		void initialize(const Statement& statement);

		/** For testing purposes */
		static OrthotropicMaterial generateRandomMaterial();
	};

}

#endif // LIBGCM_ORTHOTROPICMATERIAL_HPP
