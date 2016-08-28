#include <libgcm/rheology/models/ElasticModel.hpp>

#include <libgcm/util/math/GslUtils.hpp>

namespace gcm {

template<>
void ElasticModel<3>::
getColumnsWithRho(const int stage, int& i, int& j, int& k) {
/// indices of columns which contains (- 1 / rho) in increasing order
	switch (stage) {
	case 0 :
		i = 3; j = 4; k = 5;
		break;
	case 1 :
		i = 4; j = 6; k = 7;
		break;
	case 2 :
		i = 5; j = 7; k = 8;
		break;
	default:
		THROW_INVALID_ARG("Invalid stage number");
	}
}


template<>
void ElasticModel<3>::
getZeroColumns(const int stage, int& i, int& j, int& k) {
/// indices of columns with all zeros in increasing order
	switch (stage) {
	case 0 :
		i = 6; j = 7; k = 8;
		break;
	case 1 :
		i = 3; j = 5; k = 8;
		break;
	case 2 :
		i = 3; j = 4; k = 6;
		break;
	default:
		THROW_INVALID_ARG("Invalid stage number");
	}
}


template<>
linal::VECTOR<3, long double> ElasticModel<3>::
constructEigenvaluesPolynomial(ConstGcmMatricesPtr m, const int stage) {
/// return coefficients of polynomial \f$ x^3 + ax^2 + bx + c $\f,
/// which roots are squared eigenvalues of the matrix A
/// the order in vector is (a, b, c)

	const auto& A = (*m)(stage).A;
	linal::VECTOR<3, long double> p;
	
	int i, j, k;
	getColumnsWithRho(stage, i, j, k);
	
	const long double r = A(0,i);
	p(0) = r * (-A(k,2) - A(j,1) - A(i,0));
	
	p(1) = r * r * ((A(j,1) + A(i,0)) * A(k,2) - A(j,2) * A(k,1) - 
		A(i,2) * A(k,0) + A(i,0) * A(j,1) - A(i,1) * A(j,0));
	
	p(2) = r * r * r * ((-A(i,0) * A(j,1) + A(i,1) * A(j,0)) * A(k,2) + 
		(A(i,0) * A(j,2) - A(i,2) * A(j,0)) * A(k,1) +
		(-A(i,1) * A(j,2) + A(i,2) * A(j,1)) * A(k,0));
	
	return p;
}


template<>
std::vector<linal::VECTOR<9, long double>> ElasticModel<3>::
findEigenvectors(const long double l, const linal::Matrix<9, 9>& A,
                 const int stage, const int numberOfVectorsToSearch) {
/// @param l eigenvalue
	
	int i, j, k; //< indices of strings used at simplified matrix
	getColumnsWithRho(stage, i, j, k);
	int p, q, m; //< indices of strings filled at the end of calculations
	getZeroColumns(stage, p, q, m);
	
	const long double r = A(0, i);
	const linal::MATRIX<3, 3, long double> M = { // simplified matrix
			A(i, 0) - l*l / r,  A(i, 1),            A(i, 2),
			A(j, 0),            A(j, 1) - l*l / r,  A(j, 2),
			A(k, 0),            A(k, 1),            A(k, 2) - l*l / r,
	};
	
	auto solution3x3 = linal::solveDegenerateLinearSystem(M, numberOfVectorsToSearch);
	
	std::vector<linal::VECTOR<9, long double>> result;
	for (const auto s3x3 : solution3x3) {
		linal::VECTOR<9, long double> ev;
		ev(0) = s3x3(0); ev(1) = s3x3(1); ev(2) = s3x3(2);
		
		ev(i) = l / r * ev(0);
		ev(j) = l / r * ev(1);
		ev(k) = l / r * ev(2);
		ev(p) = (A(p,0) * ev(0) + A(p,1) * ev(1) + A(p,2) * ev(2)) / l;
		ev(q) = (A(q,0) * ev(0) + A(q,1) * ev(1) + A(q,2) * ev(2)) / l;
		ev(m) = (A(m,0) * ev(0) + A(m,1) * ev(1) + A(m,2) * ev(2)) / l;
		
		result.push_back(ev);
	}
	
	return result;
}


template<>
std::vector<linal::VECTOR<9, long double>> ElasticModel<3>::
findEigenstrings(const long double l, const linal::Matrix<9, 9>& A,
                 const int stage, const int numberOfStringsToSearch) {
/// @param l eigenvalue
		
	int i, j, k; //< indices of strings used at simplified matrix
	getColumnsWithRho(stage, i, j, k);
	int p, q, m; //< indices of result components to fill with zeros
	getZeroColumns(stage, p, q, m);
	
	const long double r = A(0, i);
	const linal::MATRIX<3, 3, long double> M = { // simplified matrix
			A(i, 0) - l*l / r,  A(j, 0),            A(k, 0),
			A(i, 1),            A(j, 1) - l*l / r,  A(k, 1),
			A(i, 2),            A(j, 2),            A(k, 2) - l*l / r,
	};
	
	auto solution3x3 = linal::solveDegenerateLinearSystem(M, numberOfStringsToSearch);
	
	std::vector<linal::VECTOR<9, long double>> result;
	for (const auto s3x3 : solution3x3) {
		linal::VECTOR<9, long double> ev;
		ev(i) = s3x3(0); ev(j) = s3x3(1); ev(k) = s3x3(2);
		
		ev(0) = l / r * ev(i);
		ev(1) = l / r * ev(j);
		ev(2) = l / r * ev(k);
		ev(p) = 0;
		ev(q) = 0;
		ev(m) = 0;
		
		result.push_back(ev);
	}
	
	return result;
}


template<>
void ElasticModel<3>::
constructRotated(
		GcmMatricesPtr m, std::shared_ptr<const OrthotropicMaterial> material) {
/// if main axes of orthotropic material are not equal 
/// to coordinate system, this is a general case of 
/// arbitrary anisotropic material matrix decomposition
	
	const real rho = material->rho;
	const auto c = material->getRotatedElasticMatrix();
	m->clear();
	
	m->m[0].A = {
		 0,        0,         0,         -1.0 / rho, 0, 0, 0, 0, 0,
		 0,        0,         0,         0, -1.0 / rho, 0, 0, 0, 0,
		 0,        0,         0,         0, 0, -1.0 / rho, 0, 0, 0,
		
		-c(0, 0), -c(0, 5), -c(0, 4),    0, 0, 0, 0, 0, 0,
		-c(0, 5), -c(5, 5), -c(4, 5),    0, 0, 0, 0, 0, 0,
		-c(0, 4), -c(4, 5), -c(4, 4),    0, 0, 0, 0, 0, 0,
		-c(0, 1), -c(1, 5), -c(1, 4),    0, 0, 0, 0, 0, 0,
		-c(0, 3), -c(3, 5), -c(3, 4),    0, 0, 0, 0, 0, 0,
		-c(0, 2), -c(2, 5), -c(2, 4),    0, 0, 0, 0, 0, 0,
	};
	
	
	m->m[1].A = {
		 0,        0,         0,         0, -1.0 / rho, 0, 0, 0, 0,
		 0,        0,         0,         0, 0, 0, -1.0 / rho, 0, 0,
		 0,        0,         0,         0, 0, 0, 0, -1.0 / rho, 0,
		
		-c(0, 5), -c(0, 1), -c(0, 3),    0, 0, 0, 0, 0, 0,
		-c(5, 5), -c(1, 5), -c(3, 5),    0, 0, 0, 0, 0, 0,
		-c(4, 5), -c(1, 4), -c(3, 4),    0, 0, 0, 0, 0, 0,
		-c(1, 5), -c(1, 1), -c(1, 3),    0, 0, 0, 0, 0, 0,
		-c(3, 5), -c(1, 3), -c(3, 3),    0, 0, 0, 0, 0, 0,
		-c(2, 5), -c(1, 2), -c(2, 3),    0, 0, 0, 0, 0, 0,
	};
	
	
	m->m[2].A = {
		 0,        0,         0,         0, 0, -1.0 / rho, 0, 0, 0,
		 0,        0,         0,         0, 0, 0, 0, -1.0 / rho, 0,
		 0,        0,         0,         0, 0, 0, 0, 0, -1.0 / rho,
		
		-c(0, 4), -c(0, 3), -c(0, 2),    0, 0, 0, 0, 0, 0,
		-c(4, 5), -c(3, 5), -c(2, 5),    0, 0, 0, 0, 0, 0,
		-c(4, 4), -c(3, 4), -c(2, 4),    0, 0, 0, 0, 0, 0,
		-c(1, 4), -c(1, 3), -c(1, 2),    0, 0, 0, 0, 0, 0,
		-c(3, 4), -c(3, 3), -c(2, 3),    0, 0, 0, 0, 0, 0,
		-c(2, 4), -c(2, 3), -c(2, 2),    0, 0, 0, 0, 0, 0,
	};
	
	
	for (int stage = 0; stage < 3; stage++) {
		auto& mat = m->m[stage];
		const auto squaredEigenvalues = GslUtils::solveThirdOrderPolynomial(
				constructEigenvaluesPolynomial(m, stage));
		
		// s1-waves, s2-waves, p-waves
		real s1w = (real) std::sqrt(squaredEigenvalues(2));
		real s2w = (real) std::sqrt(squaredEigenvalues(1));
		real pw  = (real) std::sqrt(squaredEigenvalues(0));
		mat.L = {-s1w, s1w, -s2w, s2w, -pw, pw, 0, 0, 0};
			
		if (squaredEigenvalues(1) != squaredEigenvalues(2)) {
		// three different in modulus eigenvalues (s1-, s2-, p-waves)
			for (int i = 0; i < 6; i++) {
				auto evec = findEigenvectors(mat.L(i), mat.A, stage, 1);
				mat.U1.setColumn(i, evec[0]);
				auto estr = findEigenstrings(mat.L(i), mat.A, stage, 1);
				mat.U.setRow(i, estr[0]);
			}
		
		} else {
		// two eigenvalues are equal in modulus (s-, s-, p-waves)
			for (int i = 4; i < 6; i++) {
				auto evec = findEigenvectors(mat.L(i), mat.A, stage, 1);
				mat.U1.setColumn(i, evec[0]);
				auto estr = findEigenstrings(mat.L(i), mat.A, stage, 1);
				mat.U.setRow(i, estr[0]);
			}
			for (int i = 0; i < 2; i++) {
				auto evecs = findEigenvectors(mat.L(i), mat.A, stage, 2);
				mat.U1.setColumn(i,     evecs[0]);
				mat.U1.setColumn(i + 2, evecs[1]);
				auto estr = findEigenstrings(mat.L(i), mat.A, stage, 2);
				mat.U.setRow(i,     estr[0]);
				mat.U.setRow(i + 2, estr[1]);
			}
		}
		
		// set eigenvectors for zero eigenvalues
		int p, q, r;
		getZeroColumns(stage, p, q, r);
		mat.U1(p, 6) = mat.U1(q, 7) = mat.U1(r, 8) = 1;
		
		// set eigenstrings for zero eigenvalues
		mat.U(6, p) = mat.U(7, q) = mat.U(8, r) = 1;

		int i, j, k;
		getColumnsWithRho(stage, i, j, k);
		const linal::MATRIX<3, 3, long double> M = {
				mat.A(i, 0),  mat.A(j, 0),  mat.A(k, 0),
				mat.A(i, 1),  mat.A(j, 1),  mat.A(k, 1),
				mat.A(i, 2),  mat.A(j, 2),  mat.A(k, 2),
		};
		
		Real3 x = linal::solveLinearSystem(
				M, Real3({ -mat.A(p, 0), -mat.A(p, 1), -mat.A(p, 2) }));
		mat.U(6, i) = x(0); mat.U(6, j) = x(1); mat.U(6, k) = x(2);
		
		x = linal::solveLinearSystem(
				M, Real3({ -mat.A(q, 0), -mat.A(q, 1), -mat.A(q, 2) }));
		mat.U(7, i) = x(0); mat.U(7, j) = x(1); mat.U(7, k) = x(2);
		
		x = linal::solveLinearSystem(
				M, Real3({ -mat.A(r, 0), -mat.A(r, 1), -mat.A(r, 2) }));
		mat.U(8, i) = x(0); mat.U(8, j) = x(1); mat.U(8, k) = x(2);
		
		// now, U*U1 is diagonal matrix
		// normalize U and U1 in order U*U1 to be identity matrix
		auto diag = linal::diagonalMultiply(mat.U, mat.U1);
		for (int n = 0; n < 9; n++) {
			real normalizer = std::sqrt(std::fabs(diag(n)));
			mat.U1.setColumn(n, mat.U1.getColumn(n) / normalizer);
			mat.U.setRow(n, Utils::sign(diag(n)) * mat.U.getRow(n) / normalizer);
		}
	}
	
	m->checkDecomposition(1e-2);
}


template<>
void ElasticModel<3>::
constructNotRotated(GcmMatricesPtr m, const real rho,
		const real c11, const real c12, const real c13,
		const real c22, const real c23, const real c33,
		const real c44, const real c55, const real c66) {
/// if main axes of orthotropic material are equal 
/// to coordinate system, decomposition is trivial

	m->m[0].A = {
		0, 0, 0, -1.0 / rho, 0, 0, 0, 0, 0,
		0, 0, 0, 0, -1.0 / rho, 0, 0, 0, 0,
		0, 0, 0, 0, 0, -1.0 / rho, 0, 0, 0,
		-c11, 0, 0, 0, 0, 0, 0, 0, 0,
		0, -c66, 0, 0, 0, 0, 0, 0, 0,
		0, 0, -c55, 0, 0, 0, 0, 0, 0,
		-c12, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0,
		-c13, 0, 0, 0, 0, 0, 0, 0, 0
	};

	m->m[0].L =
	{
		-sqrt(c66 / rho), sqrt(c66 / rho), -sqrt(c55 / rho), sqrt(c55 / rho),
		-sqrt(c11 / rho), sqrt(c11 / rho), 0, 0, 0
	};

	m->m[0].U = {
		0, 1.0, 0, 0, 1.0 / (sqrt(c66) * sqrt(rho)), 0, 0, 0, 0,
		0, 1.0, 0, 0, -1.0 / (sqrt(c66) * sqrt(rho)), 0, 0, 0, 0,
		0, 0, 1.0, 0, 0, 1.0 / (sqrt(c55) * sqrt(rho)), 0, 0, 0,
		0, 0, 1.0, 0, 0, -1.0 / (sqrt(c55) * sqrt(rho)), 0, 0, 0,
		1.0, 0, 0, 1.0 / (sqrt(c11) * sqrt(rho)), 0, 0, 0, 0, 0,
		1.0, 0, 0, -1.0 / (sqrt(c11) * sqrt(rho)), 0, 0, 0, 0, 0,
		0, 0, 0, -c12 / c11, 0, 0, 1.0, 0, 0.0,
		0, 0, 0, 0, 0, 0, 0, 1.0, 0,
		0, 0, 0, -(1.0 * c13) / c11, 0, 0, 0, 0, 1.0
	};

	m->m[0].U1 = {
		0, 0, 0, 0, 0.5, 0.5, 0, 0, 0,
		0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0.5, 0.5, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0.5 * sqrt(c11) * sqrt(rho), -0.5 * sqrt(c11) * sqrt(rho), 0, 0, 0,
		0.5 * sqrt(c66) * sqrt(rho), -0.5 * sqrt(c66) * sqrt(rho), 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0.5 * sqrt(c55) * sqrt(rho), -0.5 * sqrt(c55) * sqrt(rho), 0, 0, 0, 0, 0,
		0, 0, 0, 0, (0.5 * c12 * sqrt(rho)) / sqrt(c11),
		-(0.5 * c12 * sqrt(rho)) / sqrt(c11), 1, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 1, 0,
		0, 0, 0, 0, (0.5 * c13 * sqrt(rho)) / sqrt(c11),
		-(0.5 * c13 * sqrt(rho)) / sqrt(c11), 0, 0, 1
	};


	m->m[1].A = {
		0, 0, 0, 0, -1.0 / rho, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, -1.0 / rho, 0, 0,
		0, 0, 0, 0, 0, 0, 0, -1.0 / rho, 0,
		0, -c12, 0, 0, 0, 0, 0, 0, 0,
		-c66, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, -c22, 0, 0, 0, 0, 0, 0, 0,
		0, 0, -c44, 0, 0, 0, 0, 0, 0,
		0, -c23, 0, 0, 0, 0, 0, 0, 0
	};

	m->m[1].L =
	{
		-sqrt(c66 / rho), sqrt(c66 / rho), -sqrt(c44 / rho), sqrt(c44 / rho),
		-sqrt(c22 / rho), sqrt(c22 / rho), 0, 0, 0
	};

	m->m[1].U = {
		1.0, 0, 0, 0, 1.0 / (sqrt(c66) * sqrt(rho)), 0, 0, 0, 0,
		1.0, 0, 0, 0, -1.0 / (sqrt(c66) * sqrt(rho)), 0, 0, 0, 0,
		0, 0, 1.0, 0, 0, 0, 0, 1.0 / (sqrt(c44) * sqrt(rho)), 0,
		0, 0, 1.0, 0, 0, 0, 0, -1.0 / (sqrt(c44) * sqrt(rho)), 0,
		0, 1.0, 0, 0, 0, 0, 1.0 / (sqrt(c22) * sqrt(rho)), 0, 0,
		0, 1.0, 0, 0, 0, 0, -1.0 / (sqrt(c22) * sqrt(rho)), 0, 0,
		0, 0, 0, 1.0, 0, 0, -c12 / c22, 0, 0.0,
		0, 0, 0, 0, 0, 1.0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, -(1.0 * c23) / c22, 0, 1.0
	};

	m->m[1].U1 = {
		0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0.5, 0.5, 0, 0, 0,
		0, 0, 0.5, 0.5, 0, 0, 0, 0, 0,
		0, 0, 0, 0, (0.5 * c12) / sqrt(c22 / rho),
		-(0.5 * c12) / sqrt(c22 / rho), 1, 0, 0,
		0.5 * sqrt(c66) * sqrt(rho), -0.5 * sqrt(c66) * sqrt(rho), 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 1, 0,
		0, 0, 0, 0, 0.5 * sqrt(c22) * sqrt(rho), -0.5 * sqrt(c22) * sqrt(rho), 0, 0, 0,
		0, 0, 0.5 * sqrt(c44) * sqrt(rho), -0.5 * sqrt(c44) * sqrt(rho), 0, 0, 0, 0, 0,
		0, 0, 0, 0, (0.5 * c23 * sqrt(rho)) / sqrt(c22),
		-(0.5 * c23 * sqrt(rho)) / sqrt(c22), 0, 0, 1
	};


	m->m[2].A = {
		0, 0, 0, 0, 0, -1.0 / rho, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, -1.0 / rho, 0,
		0, 0, 0, 0, 0, 0, 0, 0, -1.0 / rho,
		0, 0, -c13, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0,
		-c55, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, -c23, 0, 0, 0, 0, 0, 0,
		0, -c44, 0, 0, 0, 0, 0, 0, 0,
		0, 0, -c33, 0, 0, 0, 0, 0, 0
	};

	m->m[2].L =
	{
		-sqrt(c55 / rho), sqrt(c55 / rho), -sqrt(c44 / rho), sqrt(c44 / rho),
		-sqrt(c33 / rho), sqrt(c33 / rho), 0, 0, 0
	};

	m->m[2].U = {
		1.0, 0, 0, 0, 0, 1.0 / (sqrt(c55) * sqrt(rho)), 0, 0, 0,
		1.0, 0, 0, 0, 0, -1.0 / (sqrt(c55) * sqrt(rho)), 0, 0, 0,
		0, 1.0, 0, 0, 0, 0, 0, 1.0 / (sqrt(c44) * sqrt(rho)), 0,
		0, 1.0, 0, 0, 0, 0, 0, -1.0 / (sqrt(c44) * sqrt(rho)), 0,
		0, 0, 1.0, 0, 0, 0, 0, 0, 1.0 / (sqrt(c33) * sqrt(rho)),
		0, 0, 1.0, 0, 0, 0, 0, 0, -1.0 / (sqrt(c33) * sqrt(rho)),
		0, 0, 0, 1.0, 0, 0, 0, 0, -(1.0 * c13) / c33,
		0, 0, 0, 0, 1.0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 1.0, 0, -(1.0 * c23) / c33
	};

	m->m[2].U1 = {
		0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0.5, 0.5, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0.5, 0.5, 0, 0, 0,
		0, 0, 0, 0, (0.5 * c13 * sqrt(rho)) / sqrt(c33),
		-(0.5 * c13 * sqrt(rho)) / sqrt(c33), 1, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 1, 0,
		0.5 * sqrt(c55) * sqrt(rho), -0.5 * sqrt(c55) * sqrt(rho), 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, (0.5 * c23 * sqrt(rho)) / sqrt(c33),
		-(0.5 * c23 * sqrt(rho)) / sqrt(c33), 0, 0, 1,
		0, 0, 0.5 * sqrt(c44) * sqrt(rho), -0.5 * sqrt(c44) * sqrt(rho), 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0.5 * sqrt(c33) * sqrt(rho), -0.5 * sqrt(c33) * sqrt(rho), 0, 0, 0
	};

}


template<>
void ElasticModel<3>::constructGcmMatrices(GcmMatricesPtr m,
		std::shared_ptr<const OrthotropicMaterial> material,
		const MatrixDD& basis) {
	assert_true(basis == linal::identity(basis)); // TODO
	
	if (material->anglesOfRotation == Real3::Zeros()) {
		constructNotRotated(m, material->rho, 
				material->c11, material->c12, material->c13,
				material->c22, material->c23, material->c33,
				material->c44, material->c55, material->c66);
		
		m->checkDecomposition();
	
	} else {
		constructRotated(m, material);
		
	}
	
}


}
