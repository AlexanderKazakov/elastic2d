#ifndef LIBGCM_ORTHOTROPICMATERIAL_HPP
#define LIBGCM_ORTHOTROPICMATERIAL_HPP

#include <libgcm/rheology/materials/AbstractMaterial.hpp>

namespace gcm {
struct IsotropicMaterial;

struct OrthotropicMaterial : public AbstractMaterial {
	static const Materials::T Type = Materials::T::ORTHOTROPIC;
	
	real rho = 0;          ///< density
	
	union {
		real c[9];         ///< elastic coefficients
			struct {
				real c11, c12, c13,
				          c22, c23,
				               c33,
				                    c44,
				                         c55,
				                              c66;
			};
	};
	
	real yieldStrength = 0;            ///< plasticity parameters
	real continualDamageParameter = 0; ///< parameter in continual damage equation
	real tau0 = 0;                     ///< viscosity parameter (decay time
	
	/// basis of elastic constants (main axes of material)
	/// is rotated relative to global basis by this angles
	Real3 anglesOfRotation = Real3::Zeros();

	OrthotropicMaterial(const OrthotropicMaterial& other) = default;
	OrthotropicMaterial(const IsotropicMaterial& isotropic);
	OrthotropicMaterial(const real rho_ = 0, std::initializer_list<real> =
	                    {0, 0, 0, 0, 0, 0, 0, 0, 0},
	                    const real yieldStrength_ = 0,
	                    const real continualDamageParameter_ = 0,
	                    const Real3 phi = Real3::Zeros(),
	                    const real tau0_ = 0);
	
	/** Return matrix of elastic coefficients in main axes of material */
	ElasticMatrix getElasticMatrix() const {
		auto ans = ElasticMatrix::Zeros();
		
		ans(0, 0) = c11; ans(0, 1) = c12; ans(0, 2) = c13;
		                 ans(1, 1) = c22; ans(1, 2) = c23;
		                                  ans(2, 2) = c33;
		
		ans(3, 3) = c44; ans(4, 4) = c55; ans(5, 5) = c66;
		
		return ans;
	}
	
	/** Return matrix of elastic coefficients in global coordinate system */
	ElasticMatrix getRotatedElasticMatrix() const {
		return rotate(getElasticMatrix(), anglesOfRotation);
	}
	
	/** For testing purposes */
	static OrthotropicMaterial generateRandomMaterial(const bool rotate = false);

};
}

#endif // LIBGCM_ORTHOTROPICMATERIAL_HPP
