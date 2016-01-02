#include <gtest/gtest.h>

#include "lib/config.hpp"
#include "lib/PDEMatrices.hpp"

using namespace gcm;
using namespace gcm::linal;

const int NUMBER_ITERATIONS = 1000;

const real RHO_MAX = 100.0;
const real RHO_MIN = 0.01;
const real LAMBDA_MAX = 1e+6;
const real LAMBDA_MIN = 1.0;
const real MU_MAX = 1e+6;
const real MU_MIN = 1.0;

const real rho0 = 8.0; // default density
const real lambda0 = 12e+4; // default Lame parameter
const real mu0 = 77e+3; // default Lame parameter


TEST(Linal, MatrixMatrixMultiplication)
{
	Matrix<5,5> A; Matrix<5,5> B; Matrix<5,5> C; // C = A * B
	
	A(0, 0) = 0;   A(0, 1) = 1;   A(0, 2) = 2;   A(0, 3) = 3;   A(0, 4) = 0;
	A(1, 0) = 4;   A(1, 1) = 5;   A(1, 2) = 6;   A(1, 3) = 0;   A(1, 4) = 0;
	A(2, 0) = 0;   A(2, 1) = 0;   A(2, 2) = 1;   A(2, 3) = 0;   A(2, 4) = 1;
	A(3, 0) = -2;  A(3, 1) = 5;   A(3, 2) = 4;   A(3, 3) = -7;  A(3, 4) = 0;
	A(4, 0) = 5;   A(4, 1) = 6;   A(4, 2) = 8;   A(4, 3) = 0;   A(4, 4) = 0;

	B(0, 0) = 2;   B(0, 1) = 0;   B(0, 2) = 0;   B(0, 3) = 0;   B(0, 4) = 5;
	B(1, 0) = 0;   B(1, 1) = 0;   B(1, 2) = 1;   B(1, 3) = 0;   B(1, 4) = 1;
	B(2, 0) = -6;  B(2, 1) = 0;   B(2, 2) = 0;   B(2, 3) = 0;   B(2, 4) = 0;
	B(3, 0) = 0;   B(3, 1) = -7;  B(3, 2) = 0;   B(3, 3) = -8;  B(3, 4) = 0;
	B(4, 0) = 0;   B(4, 1) = 0;   B(4, 2) = -11; B(4, 3) = 0;   B(4, 4) = 0;

	C(0, 0) = -12; C(0, 1) = -21; C(0, 2) = 1;   C(0, 3) = -24; C(0, 4) = 1;
	C(1, 0) = -28; C(1, 1) = 0;   C(1, 2) = 5;   C(1, 3) = 0;   C(1, 4) = 25;
	C(2, 0) = -6;  C(2, 1) = 0;   C(2, 2) = -11; C(2, 3) = 0;   C(2, 4) = 0;
	C(3, 0) = -28; C(3, 1) = 49;  C(3, 2) = 5;   C(3, 3) = 56;  C(3, 4) = -5;
	C(4, 0) = -38; C(4, 1) = 0;   C(4, 2) = 6;   C(4, 3) = 0;   C(4, 4) = 31;

	Matrix<5,5> AB = A * B;
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			ASSERT_NEAR(AB(i, j), C(i, j), EQUALITY_TOLERANCE);
		}
	}
}


TEST(Linal, MatrixVectorMultiplication)
{
	Matrix<5, 5> A; Vector<5> b; Vector<5> c; // c = A * b
	
	A(0, 0) = 0;   A(0, 1) = 1;   A(0, 2) = 2;   A(0, 3) = 3;   A(0, 4) = 0;
	A(1, 0) = 4;   A(1, 1) = 5;   A(1, 2) = 6;   A(1, 3) = 0;   A(1, 4) = 0;
	A(2, 0) = 0;   A(2, 1) = 0;   A(2, 2) = 1;   A(2, 3) = 0;   A(2, 4) = 1;
	A(3, 0) = -2;  A(3, 1) = 5;   A(3, 2) = 4;   A(3, 3) = -7;  A(3, 4) = 0;
	A(4, 0) = 5;   A(4, 1) = 6;   A(4, 2) = 8;   A(4, 3) = 0;   A(4, 4) = 0;
	
	b(0) = 2;      b(1) = 0;      b(2) = -13;    b(3) = 0;      b(4) = 1;

	c(0) = -26;    c(1) = -70;    c(2) = -12;    c(3) = -56;    c(4) = -94;

	Vector<5> Ab = A * b;
	for (int i = 0; i < 5; i++) {
		ASSERT_NEAR(Ab(i), c(i), EQUALITY_TOLERANCE);
	}
}


TEST(Linal, TraceVerification)
{
	Matrix A;
	A(0, 0) = 12; A(1, 1) = 56.333; A(2, 2) = 1; A(3, 3) = 0; A(4, 4) = -34.0022;
	ASSERT_NEAR(A.trace(), 35.3308, EQUALITY_TOLERANCE);
}


TEST(Linal, TraceComparison)
{
	PDEMatrices matrix(rho0, lambda0, mu0);
	ASSERT_NEAR(matrix.A(0).A.trace(), matrix.A(0).L.trace(), EQUALITY_TOLERANCE)
		<< "(x) A = " << matrix.A(0).A << "L = " << matrix.A(0).L;
	ASSERT_NEAR(matrix.A(1).A.trace(), matrix.A(1).L.trace(), EQUALITY_TOLERANCE)
		<< "(y) A = " << matrix.A(1).A << "L = " << matrix.A(1).L;

	srand(time(0));
	for(int i = 0; i < NUMBER_ITERATIONS; i++) {
		double rho = ((RHO_MAX - RHO_MIN) * rand()) / RAND_MAX + RHO_MIN;
		real lambda = ((LAMBDA_MAX - LAMBDA_MIN) * rand()) / RAND_MAX + LAMBDA_MIN;
		real mu = ((MU_MAX - MU_MIN) * rand()) / RAND_MAX + MU_MIN;
		PDEMatrices matrix1(rho, lambda, mu);
		ASSERT_NEAR(matrix1.A(0).A.trace(), matrix1.A(0).L.trace(), EQUALITY_TOLERANCE)
			<< "(x) A = " << matrix1.A(0).A << "L = " << matrix1.A(0).L;
		ASSERT_NEAR(matrix1.A(1).A.trace(), matrix1.A(1).L.trace(), EQUALITY_TOLERANCE)
			<< "(y) A = " << matrix1.A(1).A << "L = " << matrix1.A(1).L;
	}
}


TEST(Linal, LeftEigenVectors)
{
	PDEMatrices matrix(rho0, lambda0, mu0);
	Matrix AU1 = matrix.A(0).A * matrix.A(0).U1;
	Matrix U1L = matrix.A(0).U1 * matrix.A(0).L;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			ASSERT_NEAR(AU1(i, j), U1L(i, j), EQUALITY_TOLERANCE) << "(x) U1 = " << matrix.A(0).U1;
		}
	}
	AU1 = matrix.A(1).A * matrix.A(1).U1;
	U1L = matrix.A(1).U1 * matrix.A(1).L;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			ASSERT_NEAR(AU1(i, j), U1L(i, j), EQUALITY_TOLERANCE) << "(y) U1 = " << matrix.A(1).U1;
		}
	}

	srand(time(0));
	for(int i = 0; i < NUMBER_ITERATIONS; i++) {
		real rho = ((RHO_MAX - RHO_MIN) * rand()) / RAND_MAX + RHO_MIN;
		real lambda = ((LAMBDA_MAX - LAMBDA_MIN) * rand()) / RAND_MAX + LAMBDA_MIN;
		real mu = ((MU_MAX - MU_MIN) * rand()) / RAND_MAX + MU_MIN;
		PDEMatrices matrix1(rho, lambda, mu);

		AU1 = matrix1.A(0).A * matrix1.A(0).U1;
		U1L = matrix1.A(0).U1 * matrix1.A(0).L;
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++) {
				ASSERT_NEAR(AU1(i, j), U1L(i, j), EQUALITY_TOLERANCE) << "(x) U1 = " << matrix1.A(0).U1;
			}
		}
		AU1 = matrix1.A(1).A * matrix1.A(1).U1;
		U1L = matrix1.A(1).U1 * matrix1.A(1).L;
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++) {
				ASSERT_NEAR(AU1(i, j), U1L(i, j), EQUALITY_TOLERANCE) << "(y) U1 = " << matrix1.A(1).U1;
			}
		}
	}

}


TEST(Linal, RightEigenVectors)
{
	PDEMatrices matrix(rho0, lambda0, mu0);
	Matrix UA = matrix.A(0).U * matrix.A(0).A;
	Matrix LU = matrix.A(0).L * matrix.A(0).U;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			ASSERT_NEAR(UA(i, j), LU(i, j), EQUALITY_TOLERANCE) << "(x) U = " << matrix.A(0).U;
		}
	}
	UA = matrix.A(0).U * matrix.A(0).A;
	LU = matrix.A(0).L * matrix.A(0).U;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			ASSERT_NEAR(UA(i, j), LU(i, j), EQUALITY_TOLERANCE) << "(y) U = " << matrix.A(1).U;
		}
	}

	srand(time(0));
	for(int i = 0; i < NUMBER_ITERATIONS; i++) {
		real rho = ((RHO_MAX - RHO_MIN) * rand()) / RAND_MAX + RHO_MIN;
		real lambda = ((LAMBDA_MAX - LAMBDA_MIN) * rand()) / RAND_MAX + LAMBDA_MIN;
		real mu = ((MU_MAX - MU_MIN) * rand()) / RAND_MAX + MU_MIN;
		PDEMatrices matrix1(rho, lambda, mu);

		UA = matrix1.A(0).U * matrix1.A(0).A;
		LU = matrix1.A(0).L * matrix1.A(0).U;
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++) {
				ASSERT_NEAR(UA(i, j), LU(i, j), EQUALITY_TOLERANCE) << "(x) U = " << matrix1.A(0).U;
			}
		}
		UA = matrix1.A(0).U * matrix1.A(0).A;
		LU = matrix1.A(0).L * matrix1.A(0).U;
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++) {
				ASSERT_NEAR(UA(i, j), LU(i, j), EQUALITY_TOLERANCE) << "(y) U = " << matrix1.A(1).U;
			}
		}

	}
}


TEST(Linal, InverseMatrix)
{
	PDEMatrices matrix(rho0, lambda0, mu0);
	Matrix UU1 = matrix.A(0).U * matrix.A(0).U1;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			ASSERT_NEAR(UU1(i, j), (i == j), EQUALITY_TOLERANCE) << "(x) UU1 = " << UU1;
		}
	}
	UU1 = matrix.A(1).U * matrix.A(1).U1;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			ASSERT_NEAR(UU1(i, j), (i == j), EQUALITY_TOLERANCE) << "(y) UU1 = " << UU1;
		}
	}

	srand(time(0));
	for(int i = 0; i < NUMBER_ITERATIONS; i++) {
		real rho = ((RHO_MAX - RHO_MIN) * rand()) / RAND_MAX + RHO_MIN;
		real lambda = ((LAMBDA_MAX - LAMBDA_MIN) * rand()) / RAND_MAX + LAMBDA_MIN;
		real mu = ((MU_MAX - MU_MIN) * rand()) / RAND_MAX + MU_MIN;
		PDEMatrices matrix1(rho, lambda, mu);

		UU1 = matrix1.A(0).U * matrix1.A(0).U1;
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++) {
				ASSERT_NEAR(UU1(i, j), (i == j), EQUALITY_TOLERANCE) << "(x) UU1 = " << UU1;
			}
		}
		UU1 = matrix1.A(1).U * matrix1.A(1).U1;
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++) {
				ASSERT_NEAR(UU1(i, j), (i == j), EQUALITY_TOLERANCE) << "(y) UU1 = " << UU1;
			}
		}

	}
}


TEST(Linal, setColumn)
{
	Matrix<15,33> matrix;
	for (int i = 0; i < matrix.N; i++) {
		Vector<15> vector;
		for (int j = 0; j < matrix.M; j++) {
			vector(j) = i;
		}
		matrix.setColumn(i, vector);
	}

	for (int i = 0; i < matrix.M; i++) {
		for (int j = 0; j < matrix.N; j++) {
			ASSERT_EQ(matrix(i, j), j);
		}
	}
}


TEST(Linal, getColumn)
{
	Matrix<22,13> matrix;
	for (int i = 0; i < matrix.M; i++) {
		for (int j = 0; j < matrix.N; j++) {
			matrix(i, j) = j;
		}
	}

	for (int k = 0; k < matrix.N; k++) {
		Vector<22> column = matrix.getColumn(k);
		for (int i = 0; i < column.M; i++) {
			ASSERT_EQ(column(i), k);
		}
	}
}


TEST(Linal, diagonalMultiply)
{
	Matrix<9, 9> A, B;
	srand(time(0));

	for (int l = 0; l < NUMBER_ITERATIONS; l++) {

		for (int i = 0; i < A.M; i++) {
			for (int j = 0; j < A.N; j++) {
				A(i, j) = rand() - RAND_MAX / 2.0;
				B(i, j) = rand() - RAND_MAX / 2.0;
			}
		}

		Matrix<9,9> C = A * B;
		Vector<9> d = A.diagonalMultiply(B);

		for (int k = 0; k < d.M; k++) {
			ASSERT_EQ(C(k, k), d(k));
		}

	}
}


TEST(Linal, getDiagonalMultipliedBy)
{
	Matrix<11,11> matrix;
	for (int i = 0; i < matrix.M; i++) {
		for (int j = 0; j < matrix.N; j++) {
			matrix(i, j) = (real) i;
		}
	}
	Vector<11> vector = matrix.getDiagonalMultipliedBy(- 5.0);
	for (int k = 0; k < vector.M; k++) {
		ASSERT_EQ(vector(k), - 5.0 * k);
	}
}


TEST(Linal, VectorOperators)
{
	Vector<7> vector;
	for (int i = 0; i < vector.M; i++) {
		vector(i) = (real) i - 3;
	}

	vector = vector * 2;
	for (int j = 0; j < vector.M; j++) {
		ASSERT_EQ(vector(j), 2 * ((real) j - 3));
	}

	vector += vector;
	for (int j = 0; j < vector.M; j++) {
		ASSERT_EQ(vector(j), 4 * ((real) j - 3));
	}

	vector = vector - vector * 0.5;
	for (int j = 0; j < vector.M; j++) {
		ASSERT_EQ(vector(j), 2 * ((real) j - 3));
	}
}