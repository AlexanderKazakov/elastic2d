#ifndef LIBGCM_ORTHOTROPICMATERIAL_HPP
#define LIBGCM_ORTHOTROPICMATERIAL_HPP


#include <lib/util/Types.hpp>
#include <initializer_list>

namespace gcm {
	class IsotropicMaterial;
	template<int TM, int Dimensionality, class TMaterial> class GcmMatrices;


	/**
	 * Orthotropic elastic materials
	 */
	class OrthotropicMaterial {
	public:
		real rho;
		union {
			real c[9];
			struct {
				real c11, c12, c13,
				          c22, c23,
				               c33,
				                    c44,
				                         c55,
				                              c66;
			};
		};

		OrthotropicMaterial(const real rho, std::initializer_list<real>);
		OrthotropicMaterial(const IsotropicMaterial& isotropic);

		/** Fill in gcm matrices */
		void constructGcmMatrices(GcmMatrices<9, 3, OrthotropicMaterial>& m) const;

		/** For testing purposes */
		static OrthotropicMaterial generateRandomMaterial();
	};

};

#endif // LIBGCM_ORTHOTROPICMATERIAL_HPP
