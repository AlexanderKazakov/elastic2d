#include <gtest/gtest.h>

#include "lib/config.hpp"
#include "lib/gcm_matrices/IdealElastic2DGcmMatrices.hpp"

using namespace gcm;
using namespace gcm::linal;

const int NUMBER_ITERATIONS = 1000;

const real RHO_MAX = 100.0;
const real RHO_MIN = 0.01;
const real LAMBDA_MAX = 1e+6;
const real LAMBDA_MIN = 1.0;
const real MU_MAX = 1e+6;
const real MU_MIN = 1.0;

const real rho0 = 8.0;
const real lambda0 = 12e+4;
const real mu0 = 77e+3;


template<int N, int Dimensionality>
void traceComparison(const GcmMatrices<N, Dimensionality>& matrix) {
	for (int i = 0; i < matrix.DIMENSIONALITY; i++) {
		ASSERT_NEAR(matrix.A(i).A.trace(), matrix.A(i).L.trace(), EQUALITY_TOLERANCE)
									<< "(" << i << ") A = " << matrix.A(i).A << "L = " << matrix.A(i).L;
	}
}

TEST(GcmMatrices, TraceComparison)
{
	IdealElastic2DGcmMatrices matrix(rho0, lambda0, mu0);
	traceComparison(matrix);

	srand(time(0));
	for(int i = 0; i < NUMBER_ITERATIONS; i++) {
		real rho = ((RHO_MAX - RHO_MIN) * rand()) / RAND_MAX + RHO_MIN;
		real lambda = ((LAMBDA_MAX - LAMBDA_MIN) * rand()) / RAND_MAX + LAMBDA_MIN;
		real mu = ((MU_MAX - MU_MIN) * rand()) / RAND_MAX + MU_MIN;
		matrix = IdealElastic2DGcmMatrices(rho, lambda, mu);
		traceComparison(matrix);
	}
}

template<int N, int Dimensionality>
void leftEigenVectorCheck(const GcmMatrices<N, Dimensionality>& matrix) {
	for (int s = 0; s < matrix.DIMENSIONALITY; s++) {
		Matrix<N,N> AU1 = matrix.A(s).A * matrix.A(s).U1;
		Matrix<N,N> U1L = matrix.A(s).U1 * matrix.A(s).L;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				ASSERT_NEAR(AU1(i, j), U1L(i, j), EQUALITY_TOLERANCE) << "(" << s << ") U1 = " << matrix.A(s).U1;
			}
		}
	}
}

TEST(GcmMatrices, LeftEigenVectors)
{
	IdealElastic2DGcmMatrices matrix(rho0, lambda0, mu0);
	leftEigenVectorCheck(matrix);

	srand(time(0));
	for(int i = 0; i < NUMBER_ITERATIONS; i++) {
		real rho = ((RHO_MAX - RHO_MIN) * rand()) / RAND_MAX + RHO_MIN;
		real lambda = ((LAMBDA_MAX - LAMBDA_MIN) * rand()) / RAND_MAX + LAMBDA_MIN;
		real mu = ((MU_MAX - MU_MIN) * rand()) / RAND_MAX + MU_MIN;
		matrix = IdealElastic2DGcmMatrices(rho, lambda, mu);
		leftEigenVectorCheck(matrix);
	}
}

template<int N, int Dimensionality>
void rightEigenVectorCheck(const GcmMatrices<N, Dimensionality>& matrix) {
	for (int s = 0; s < matrix.DIMENSIONALITY; s++) {
		Matrix<N,N> UA = matrix.A(s).U * matrix.A(s).A;
		Matrix<N,N> LU = matrix.A(s).L * matrix.A(s).U;
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++) {
				ASSERT_NEAR(UA(i, j), LU(i, j), EQUALITY_TOLERANCE) << "(" << s << ") U = " << matrix.A(s).U;
			}
		}
	}
}

TEST(GcmMatrices, RightEigenVectors)
{
	IdealElastic2DGcmMatrices matrix(rho0, lambda0, mu0);
	rightEigenVectorCheck(matrix);

	srand(time(0));
	for(int i = 0; i < NUMBER_ITERATIONS; i++) {
		real rho = ((RHO_MAX - RHO_MIN) * rand()) / RAND_MAX + RHO_MIN;
		real lambda = ((LAMBDA_MAX - LAMBDA_MIN) * rand()) / RAND_MAX + LAMBDA_MIN;
		real mu = ((MU_MAX - MU_MIN) * rand()) / RAND_MAX + MU_MIN;
		matrix = IdealElastic2DGcmMatrices(rho, lambda, mu);
		rightEigenVectorCheck(matrix);
	}
}

template<int N, int Dimensionality>
void inverseCheck(const GcmMatrices<N, Dimensionality>& matrix) {
	for (int s = 0; s < matrix.DIMENSIONALITY; s++) {
		Matrix<N,N> UU1 = matrix.A(s).U * matrix.A(s).U1;
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++) {
				ASSERT_NEAR(UU1(i, j), (i == j), EQUALITY_TOLERANCE) << "(" << s << ") UU1 = " << UU1;
			}
		}
	}
}

TEST(GcmMatrices, InverseMatrix)
{
	IdealElastic2DGcmMatrices matrix(rho0, lambda0, mu0);
	inverseCheck(matrix);

	srand(time(0));
	for(int i = 0; i < NUMBER_ITERATIONS; i++) {
		real rho = ((RHO_MAX - RHO_MIN) * rand()) / RAND_MAX + RHO_MIN;
		real lambda = ((LAMBDA_MAX - LAMBDA_MIN) * rand()) / RAND_MAX + LAMBDA_MIN;
		real mu = ((MU_MAX - MU_MIN) * rand()) / RAND_MAX + MU_MIN;
		matrix = IdealElastic2DGcmMatrices(rho, lambda, mu);
		inverseCheck(matrix);
	}
}
