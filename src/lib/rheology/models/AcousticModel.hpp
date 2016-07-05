#ifndef LIBGCM_ACOUSTICMODEL_HPP
#define LIBGCM_ACOUSTICMODEL_HPP

#include <lib/rheology/models/Model.hpp>


namespace gcm {


/**
 * The model where tension tensor reduced to single pressure scalar.
 */
template<int Dimensionality>
class AcousticModel {
public:
	static const int DIMENSIONALITY = Dimensionality;

	typedef AcousticVariables<DIMENSIONALITY> PdeVariables;
	typedef DummyOde                          InternalOde;
	typedef DummyCorrector                    Corrector;
	typedef typename PdeVariables::PdeVector  PdeVector;
	
	typedef linal::Vector<DIMENSIONALITY>     RealD;
	typedef linal::Matrix<DIMENSIONALITY,
	                      DIMENSIONALITY>     MatrixDD;

	static const int PDE_SIZE = PdeVector::M;

	typedef GcmMatrices<PDE_SIZE, DIMENSIONALITY> GCM_MATRICES;
	typedef typename GCM_MATRICES::GcmMatrix      GcmMatrix;
	typedef typename InternalOde::Variables       OdeVariables;
	typedef std::shared_ptr<GCM_MATRICES>         GcmMatricesPtr;
	typedef std::shared_ptr<const GCM_MATRICES>   ConstGcmMatricesPtr;

	static const MaterialsWavesMap MATERIALS_WAVES_MAP;
	
	/// Number of characteristics with slopes of the same sign.
	/// It is equal to number of outer characteristics in border node.
	static const int OUTER_NUMBER = 1;
	/// Matrix of linear border condition "B * \vec{u} = b"
	/// @see BorderCondition
	typedef linal::Matrix<OUTER_NUMBER, PDE_SIZE> BorderMatrix;
	
	
	/**
	 * Construct gcm matrices for calculation in global basis axes
	 */
	template<typename ... Args>
	static void constructGcmMatrices(GcmMatricesPtr m,
			std::shared_ptr<const IsotropicMaterial> material, const Args& ...) {
		
		for (int i = 0; i < DIMENSIONALITY; i++) {
			RealD n = RealD::Zeros();
			n(i) = 1;
			constructGcmMatrix((*m)(i), material, linal::createLocalBasis(n));
		}
	}
	
	
	/**
	 * Construct gcm matrix along the given direction.
	 * @see ElasticModel::constructGcmMatrix for details
	 */
	inline static void 
	constructGcmMatrix(
			GcmMatrix& m, std::shared_ptr<const IsotropicMaterial> material,
			const MatrixDD& basis, const real l = 1);
	
	
	/**
	 * Matrix of linear border condition in local (connected with border) 
	 * basis for the case of fixed pressure on border.
	 * @see BorderCondition
	 */
	static BorderMatrix borderMatrixFixedForce
	(const linal::Vector<DIMENSIONALITY>&) {
		PdeVariables pde = PdeVariables::Zeros();
		pde.pressure() = 1;
		
		BorderMatrix B_;
		B_.setRow(0, pde);
		
		return B_;
	}
	
	
	/**
	 * Matrix of linear border condition in local (connected with border) 
	 * basis for the case of fixed normal velocity on border.
	 * @see BorderCondition
	 */
	static BorderMatrix borderMatrixFixedVelocity
	(const linal::Vector<DIMENSIONALITY>& borderNormal) {
		PdeVariables pde = PdeVariables::Zeros();
		pde.setVelocity(borderNormal);
		
		BorderMatrix B_;
		B_.setRow(0, pde);
		
		return B_;
	}
	
	
};


template<int Dimensionality>
inline void AcousticModel<Dimensionality>:: 
constructGcmMatrix(
		GcmMatrix& m, std::shared_ptr<const IsotropicMaterial> material,
		const MatrixDD& basis, const real l) {
	
	const real rho = material->rho;
	const real lambda = material->lambda;
	const real c1 = sqrt(lambda / rho);
	
	
	RealD n[DIMENSIONALITY];
	for (int i = 0; i < DIMENSIONALITY; i++) {
		n[i] = basis.getColumn((i + DIMENSIONALITY - 1) % DIMENSIONALITY);
	}
	
	
	/// fill matrix A along direction n[0] with scale l
	PdeVariables vec;
	linal::clear(m.A);
	
	vec.pressure() = 0;
	vec.setVelocity(l * lambda * n[0]);
	m.A.setRow(DIMENSIONALITY, vec);
	vec.setVelocity(l / rho * n[0]);
	m.A.setColumn(DIMENSIONALITY, vec);
	
	
	/// fill L with eigenvalues
	linal::clear(m.L);
	m.L(0) =  l*c1;
	m.L(1) = -l*c1;
	
	
	/// fill U1 with eigenvectors
	/// p-waves
	vec.setVelocity(n[0]);
	vec.pressure() = c1 * rho;
	m.U1.setColumn(0, vec);
	vec.pressure() = -vec.pressure();
	m.U1.setColumn(1, vec);
	
	/// zero eigenvalues
	vec.pressure() = 0;
	for (int i = 1; i < DIMENSIONALITY; i++) {
		vec.setVelocity(n[i]);
		m.U1.setColumn(i + 1, vec);
	}
	
	
	/// fill U with eigenstrings
	const real alpha = 0.5; //< normalizer for U1*U == I		
	/// p-waves
	vec.setVelocity(alpha * n[0]);
	vec.pressure() = alpha / (c1 * rho);
	m.U.setRow(0, vec);
	vec.pressure() = -vec.pressure();
	m.U.setRow(1, vec);
	
	/// zero eigenvalues
	vec.pressure() = 0;
	for (int i = 1; i < DIMENSIONALITY; i++) {
		vec.setVelocity(n[i]);
		m.U.setRow(i + 1, vec);
	}
	
	
	m.checkDecomposition(100*EQUALITY_TOLERANCE);
}


}

#endif // LIBGCM_ACOUSTICMODEL_HPP
