#include <lib/rheology/models/Model.hpp>

using namespace gcm;


void Elastic3DModel::
constructGcmMatrix(GcmMatrix& m, std::shared_ptr<const IsotropicMaterial> material,
		const linal::Matrix33& basis, const real l) {
	
	const real rho = material->rho;
	const real lambda = material->lambda;
	const real mu = material->mu;
	
	const Real3 n = basis.getColumn(2);
	const Real3 n1 = basis.getColumn(0);
	const Real3 n2 = basis.getColumn(1);
	

	/// fill matrix A along direction n with scale l
	PdeVariables vec;
	linal::clear(m.A);
	
	/// 3 first strings
	for (int i = 0; i < DIMENSIONALITY; i++) {
		linal::clear(vec);
		for (int j = 0; j < DIMENSIONALITY; j++) {
			vec.sigma(i, j) = -l * n(j) / rho;
		}
		m.A.setRow(i, vec);
	}
	
	/// 3 first columns
	for (int i = 0; i < DIMENSIONALITY; i++) {
		linal::clear(vec);
		for (int j = 0; j < DIMENSIONALITY; j++) {
			vec.sigma(i, j)  = -l * mu * n(j);
		}
		for (int j = 0; j < DIMENSIONALITY; j++) {
			vec.sigma(j, j) += -l * (lambda + (i == j) * mu) * n(i);
		}
		m.A.setColumn(i, vec);
	}
	
	
	/// fill L with eigenvalues
	const real c1 = sqrt((lambda + 2*mu) / rho);
	const real c2 = sqrt(mu / rho);
	m.L = { l*c1, -l*c1, l*c2, -l*c2, l*c2, -l*c2, 0, 0, 0 };
	
	
	/// fill U1 with eigenvectors
	typedef linal::SymmetricMatrix<DIMENSIONALITY> SigmaD;
	const SigmaD I = SigmaD::Identity();
	const SigmaD N00 = linal::symmDirectProduct(n, n);
	const SigmaD N01 = linal::symmDirectProduct(n, n1);
	const SigmaD N02 = linal::symmDirectProduct(n, n2);
	const SigmaD N11 = linal::symmDirectProduct(n1, n1);
	const SigmaD N12 = linal::symmDirectProduct(n1, n2);
	const SigmaD N22 = linal::symmDirectProduct(n2, n2);
	const real alpha = 0.5; //< normalizator for U*U1 = I
	
	/// p-waves
	vec.setVelocity(alpha * n);
	vec.setSigma(-alpha / c1 * (lambda * I + 2 * mu * N00));
	m.U1.setColumn(0, vec);
	vec.setSigma(-vec.getSigma());
	m.U1.setColumn(1, vec);

	/// s-waves
	vec.setVelocity(alpha * n1);
	vec.setSigma(-2 * alpha * mu / c2 * N01);
	m.U1.setColumn(2, vec);
	vec.setSigma(-vec.getSigma());
	m.U1.setColumn(3, vec);
	
	vec.setVelocity(alpha * n2);
	vec.setSigma(-2 * alpha * mu / c2 * N02);
	m.U1.setColumn(4, vec);
	vec.setSigma(-vec.getSigma());
	m.U1.setColumn(5, vec);
	
	/// zero eigenvalues
	vec.setVelocity(Real3::Zeros());
	vec.setSigma(2 * N12);
	m.U1.setColumn(6, vec);
	vec.setSigma((N11 - N22) / 2);
	m.U1.setColumn(7, vec);
	vec.setSigma((N11 + N22) / 2);
	m.U1.setColumn(8, vec);
	
	
	/// fill U with eigenstrings
	
	/// p-waves
	vec.setVelocity(n);
	vec.setSigma((2 * N00 - linal::Diag(N00)) / (-c1 * rho));
	m.U.setRow(0, vec);
	vec.setSigma(-vec.getSigma());
	m.U.setRow(1, vec);
	
	/// s-waves
	vec.setVelocity(n1);
	vec.setSigma((2 * N01 - linal::Diag(N01)) / (-c2 * rho));
	m.U.setRow(2, vec);
	vec.setSigma(-vec.getSigma());
	m.U.setRow(3, vec);
	
	vec.setVelocity(n2);
	vec.setSigma((2 * N02 - linal::Diag(N02)) / (-c2 * rho));
	m.U.setRow(4, vec);
	vec.setSigma(-vec.getSigma());
	m.U.setRow(5, vec);
	
	/// zero eigenvalues
	vec.setVelocity(Real3::Zeros());
	vec.setSigma((2 * N12 - linal::Diag(N12)));
	m.U.setRow(6, vec);
	vec.setSigma(2 * (N11 - N22) - linal::Diag((N11 - N22)));
	m.U.setRow(7, vec);
	vec.setSigma(2 * (N11 + N22) - linal::Diag((N11 + N22)) -
			2 * lambda / (lambda + 2 * mu) * (2 * N00 - linal::Diag(N00)));
	m.U.setRow(8, vec);
	
	
	m.checkDecomposition(100*EQUALITY_TOLERANCE);
}
