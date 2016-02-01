#ifndef LIBGCM_ORTHOTROPICMATERIAL_HPP
#define LIBGCM_ORTHOTROPICMATERIAL_HPP


#include <lib/util/Types.hpp>
#include <lib/rheology/variables/VelocitySigmaVariables.hpp>
#include <initializer_list>

namespace gcm {
	class IsotropicMaterial;
	template<typename TVariables, class TMaterial> class GcmMatrices;

	class OrthotropicMaterial {
	public:
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

		OrthotropicMaterial();
		OrthotropicMaterial(const IsotropicMaterial& isotropic);
		OrthotropicMaterial(const real _rho, std::initializer_list<real>);
		OrthotropicMaterial(const real _rho, const real _yieldStrength, std::initializer_list<real>);
		OrthotropicMaterial(const real _rho, const real _yieldStrength, const real _continualDamageParameter,
		                    std::initializer_list<real>);

		/** Fill in gcm matrices */
		void constructGcmMatrices(GcmMatrices<VelocitySigmaVariables<3>, OrthotropicMaterial>& m) const;

		/** For testing purposes */
		static OrthotropicMaterial generateRandomMaterial();
	};

};

#endif // LIBGCM_ORTHOTROPICMATERIAL_HPP
