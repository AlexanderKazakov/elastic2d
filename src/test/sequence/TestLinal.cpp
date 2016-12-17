#include <libgcm/linal/linal.hpp>
#include <libgcm/util/math/GslUtils.hpp>

#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

const int LINAL_TEST_NUMBER_ITERATIONS = 1000;

using namespace gcm;
using namespace gcm::linal;


TEST(Linal, segmentContains) {
	Real3 p = {0, 0, 0}, q = {1, 0, 0}, r = {1, 1, 1}, m = {2, 2, 2};
	ASSERT_TRUE (segmentContains(q, p, q, 0, 0));
	ASSERT_TRUE (segmentContains(q, p, p, 0, 0));
	ASSERT_TRUE (segmentContains(p, m, r, 0, 0));
	ASSERT_FALSE(segmentContains(p, r, m, EQUALITY_TOLERANCE, EQUALITY_TOLERANCE));
	ASSERT_FALSE(segmentContains(p, m, q, EQUALITY_TOLERANCE, EQUALITY_TOLERANCE));
	
	Real2 a = {0, 0}, b = {0, 1}, c = {4, 4}, d = {3, 3};
	ASSERT_TRUE (segmentContains(a, d,  a, 0, 0));
	ASSERT_TRUE (segmentContains(a, d,  d, 0, 0));
	ASSERT_TRUE (segmentContains(a, d,  (a + d)/2, EQUALITY_TOLERANCE, EQUALITY_TOLERANCE));
	ASSERT_FALSE(segmentContains(a, d,  c, EQUALITY_TOLERANCE, EQUALITY_TOLERANCE));
	ASSERT_FALSE(segmentContains(a, d, -c, EQUALITY_TOLERANCE, EQUALITY_TOLERANCE));
	ASSERT_FALSE(segmentContains(a, d,  b, EQUALITY_TOLERANCE, EQUALITY_TOLERANCE));
	
	ASSERT_TRUE (segmentContains(Real1({5}), Real1({6}), Real1({5}), 0));
	ASSERT_FALSE(segmentContains(Real1({5}), Real1({6}), Real1({7}), EQUALITY_TOLERANCE));
	ASSERT_FALSE(segmentContains(Real1({5}), Real1({6}), Real1(Real1({4})), EQUALITY_TOLERANCE));
}


TEST(Linal, containsTriangleIn3D) {
	Real3 b = {0, 0, 5}, c = {0, 5, 0}, d = {5, 0, 0}, e = {-1, 6, 0}, m = {0, 6, -1};
	ASSERT_TRUE (triangleContains(d, b, c, d, 0, 0));
	ASSERT_TRUE (triangleContains(d, b, c, c, 0, 0));
	ASSERT_TRUE (triangleContains(d, b, c, (d + b) / 2, EQUALITY_TOLERANCE, EQUALITY_TOLERANCE));
	ASSERT_TRUE (triangleContains(d, b, c, d/2 + b/4 + c/4, EQUALITY_TOLERANCE, EQUALITY_TOLERANCE));
	ASSERT_FALSE(triangleContains(d, b, c, e, EQUALITY_TOLERANCE, EQUALITY_TOLERANCE));
	ASSERT_FALSE(triangleContains(d, b, c, m, EQUALITY_TOLERANCE, EQUALITY_TOLERANCE));
	ASSERT_FALSE(triangleContains(d, b, c, Real3({1, 1, 1}), EQUALITY_TOLERANCE, EQUALITY_TOLERANCE));
}


TEST(Linal, contains2d) {
	Real2 a = {0, 0}, b = {0, 3}, c = {3, 0}, d = {3, 3};
	Real2 q = {1, 1}, p = {2, 2}, r = {0, -1};
	ASSERT_TRUE (triangleContains(a, b, c, q, 0));
	ASSERT_FALSE(triangleContains(a, b, c,-q, 0));
	ASSERT_TRUE (triangleContains(d, b, c, p, 0));
	ASSERT_FALSE(triangleContains(a, b, c, p, 0));
	ASSERT_FALSE(triangleContains(d, b, c, q, 0));
	ASSERT_TRUE (triangleContains(a, b, c, a, EQUALITY_TOLERANCE));
	ASSERT_TRUE (triangleContains(a, b, d, p, EQUALITY_TOLERANCE));
	ASSERT_TRUE (triangleContains(a, b, d, q, EQUALITY_TOLERANCE));
	
	ASSERT_TRUE (angleContains(a, b, c, q, 0));
	ASSERT_FALSE(angleContains(a, b, c,-q, 0));
	ASSERT_FALSE(angleContains(a, b, c, r, 0));
	ASSERT_TRUE (angleContains(a, b, c, p, 0));
	ASSERT_TRUE (angleContains(b, a, c, q, 0));
	ASSERT_FALSE(angleContains(b, a, c, p, 0));
	ASSERT_TRUE (angleContains(b, d, c, p, 0));
	ASSERT_FALSE(angleContains(b, d, c, q, 0));
	ASSERT_TRUE (angleContains(b, a, c, b, EQUALITY_TOLERANCE));
	ASSERT_TRUE (angleContains(b, a, c, a, EQUALITY_TOLERANCE));
	ASSERT_TRUE (angleContains(a, b, d, q, EQUALITY_TOLERANCE));
}


TEST(Linal, contains3d) {
	Real3 a = {0, 0, 0}, b = {0, 0, 5}, c = {0, 5, 0}, d = {5, 0, 0}, e = {5, 5, 5};
	Real3 q = {1, 1, 1}, p = {4, 4, 4}, r = {0, 0, -1}, s = {-1, 0, -1};
	ASSERT_TRUE (tetrahedronContains(a, b, c, d, q, 0));
	ASSERT_FALSE(tetrahedronContains(a, b, c, d,-q, 0));
	ASSERT_TRUE (tetrahedronContains(e, b, c, d, p, 0));
	ASSERT_FALSE(tetrahedronContains(a, b, c, d, p, 0));
	ASSERT_FALSE(tetrahedronContains(e, b, c, d, q, 0));
	ASSERT_TRUE (tetrahedronContains(a, b, c, d, a, EQUALITY_TOLERANCE));
	ASSERT_TRUE (tetrahedronContains(a, b, e, d, p, EQUALITY_TOLERANCE));
	ASSERT_TRUE (tetrahedronContains(a, b, e, d, q, EQUALITY_TOLERANCE));
	ASSERT_TRUE (tetrahedronContains(a, b, e, d, (a + b + d) / 3, EQUALITY_TOLERANCE));
	
	ASSERT_TRUE (solidAngleContains(a, b, c, d, q, 0));
	ASSERT_FALSE(solidAngleContains(a, b, c, d,-q, 0));
	ASSERT_FALSE(solidAngleContains(a, b, c, d, r, 0));
	ASSERT_FALSE(solidAngleContains(a, b, c, d, s, 0));
	ASSERT_TRUE (solidAngleContains(a, b, c, d, p, 0));
	ASSERT_TRUE (solidAngleContains(b, a, c, d, q, 0));
	ASSERT_FALSE(solidAngleContains(b, a, c, d, p, 0));
	ASSERT_TRUE (solidAngleContains(b, e, c, d, p, 0));
	ASSERT_FALSE(solidAngleContains(b, e, c, d, q, 0));
	ASSERT_TRUE (solidAngleContains(b, a, c, d, b, EQUALITY_TOLERANCE));
	ASSERT_TRUE (solidAngleContains(b, a, c, d, a, EQUALITY_TOLERANCE));
	ASSERT_TRUE (solidAngleContains(a, b, e, d, q, EQUALITY_TOLERANCE));
	ASSERT_TRUE (solidAngleContains(a, b, e, d, (a + b + d) / 3, EQUALITY_TOLERANCE));
}


TEST(Linal, reflectionDirection) {
	ASSERT_TRUE(linal::approximatelyEqual(
			Real2({0, -1}),
			reflectionDirection(normalize(Real2({1, 1})), Real2({1, 0}))));
	ASSERT_TRUE(linal::approximatelyEqual(
			Real2({-1, 0}),
			reflectionDirection(normalize(Real2({1, 1})), Real2({0, 1}))));
	ASSERT_TRUE(linal::approximatelyEqual(
			Real3({0, -1, 0}),
			reflectionDirection(normalize(Real3({0, 1, 1})), Real3({0, 0, 1}))));
	ASSERT_TRUE(linal::approximatelyEqual(
			Real3({0, 0, 1}),
			reflectionDirection(normalize(Real3({0, 0, 1})), Real3({0, 0, -1}))));
	ASSERT_TRUE(linal::approximatelyEqual(
			Real3({1, 0, 0}),
			reflectionDirection(normalize(Real3({0, 0, 1})), Real3({1, 0, 0}))));
	ASSERT_TRUE(linal::approximatelyEqual(
			normalize(Real3({1, 1, -1})),
			reflectionDirection(normalize(Real3({0, 0, 1})),
			                    normalize(Real3({1, 1, 1})))));
}


TEST(Linal, normMax) {
	ASSERT_EQ(7, normMax(Matrix22({
			1, 2,
			3, 4 })));
	ASSERT_EQ(11, normMax(Matrix22({
			5, 6,
			3, 4 })));
	ASSERT_EQ(19, normMax(Matrix33({
			1, 2, 3,
			3, 4, 10,
			11, 0, 8 })));
}


TEST(Linal, conditionNumber) {
	ASSERT_EQ(1, conditionNumber(Matrix33::Identity()));
	ASSERT_EQ(199*199, conditionNumber(Matrix22({
			100, 99,
			 99, 98 })));
}


TEST(Linal, invert) {
	Utils::seedRand();
	for (int i = 0; i < LINAL_TEST_NUMBER_ITERATIONS; i++) {
		
		auto m1 = random<Matrix11>(-1, 1);
		ASSERT_TRUE(approximatelyEqual(identity(m1), m1 * invert(m1)))
				<< "matrix = " << m1 << "identity = " << m1 * invert(m1);
		
		auto m2 = random<Matrix22>(-1, 1);
		ASSERT_TRUE(approximatelyEqual(identity(m2), m2 * invert(m2), 1e-5))
				<< "matrix = " << m2 << "identity = " << m2 * invert(m2);
		
		auto m3 = random<Matrix33>(-1, 1);
		ASSERT_TRUE(approximatelyEqual(identity(m3), m3 * invert(m3), 1e-5))
				<< "matrix = " << m3 << "identity = " << m3 * invert(m3);
	}
}


TEST(Linal, Assignment) {
	Matrix<1, 4> r = {1, 0, -3.4, 1.5};
	
	MatrixInt<1, 4> i = r;
	ASSERT_EQ((MatrixInt<1, 4>({1, 0, -3, 1})), i);
	
	Matrix<1, 4> r2 = i;
	ASSERT_EQ((MatrixInt<1, 4>({1, 0, -3, 1})), r2);
	
	Matrix<1, 4> r3 = r;
	ASSERT_EQ((Matrix<1, 4>({1, 0, -3.4, 1.5})), r3);
	
	ASSERT_EQ(1,    r3(0));
	ASSERT_EQ(0,    r3(1));
	ASSERT_EQ(-3.4, r3(2));
	ASSERT_EQ(1.5,  r3(3));
}


TEST(Linal, zeroOnes) {
	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 12; j++) {
			ASSERT_EQ(1, (Matrix<12, 12>::Ones()(i, j)));
			ASSERT_EQ(0, (Matrix<12, 12>::Zeros()(i, j)));
			ASSERT_EQ((i == j), (Matrix<12, 12>::Identity()(i, j)));
			
			ASSERT_EQ((i == j), (DiagonalMatrix<12>::Ones()(i, j)));
			ASSERT_EQ(0, (DiagonalMatrix<12>::Zeros()(i, j)));
			ASSERT_EQ((i == j), (DiagonalMatrix<12>::Identity()(i, j)));
			
			ASSERT_EQ(1, (MatrixInt<12, 12>::Ones()(i, j)));
			ASSERT_EQ(0, (MatrixInt<12, 12>::Zeros()(i, j)));
			ASSERT_EQ((i == j), (MatrixInt<12, 12>::Identity()(i, j)));
		}
	}
	
	Matrix<2, 4> m = {1, 1, 2, 3, 3, 5, 6, 0};
	auto m2 = zeros(m);
	ASSERT_EQ((Matrix<2, 4>({0, 0, 0, 0, 0, 0, 0, 0})), m2);
	
	auto m3 = ones(m);
	ASSERT_EQ((Matrix<2, 4>({1, 1, 1, 1, 1, 1, 1, 1})), m3);
	
	auto m4 = clear(m3);
	ASSERT_EQ(m2, m4);
	ASSERT_EQ(m2, m3);
	
	Matrix<2, 2> identityMatrix = identity(Matrix22());
	ASSERT_EQ((Matrix<2, 2>({1, 0, 0, 1})), identityMatrix);
}


TEST(Linal, getIndex) {
	ASSERT_EQ((SymmProps<Symmetric, 3, 3>::getIndex(0, 0)), 0);
	ASSERT_EQ((SymmProps<Symmetric, 3, 3>::getIndex(0, 2)), 2);
	ASSERT_EQ((SymmProps<Symmetric, 3, 3>::getIndex(2, 2)), 5);
	ASSERT_EQ((SymmProps<Symmetric, 5, 5>::getIndex(4, 3)), 13);

	ASSERT_EQ((SymmProps<NonSymmetric, 3, 3>::getIndex(0, 0)), 0);
	ASSERT_EQ((SymmProps<NonSymmetric, 3, 3>::getIndex(0, 2)), 2);
	ASSERT_EQ((SymmProps<NonSymmetric, 3, 3>::getIndex(2, 2)), 8);
	ASSERT_EQ((SymmProps<NonSymmetric, 5, 5>::getIndex(4, 3)), 23);

	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			ASSERT_EQ((SymmProps<Symmetric, 10, 10>::getIndex(i, j)),
			          (SymmProps<Symmetric, 10, 10>::getIndex(j, i)));
		}
	}
}


TEST(Linal, MatrixConstruct) {

	Matrix33 m2({1.0, 2.0, 3.0,
	             4.0, 5.0, 6.0,
	             7.0, 8.0, 9.0});

	ASSERT_EQ(m2(0, 0), 1.0);
	ASSERT_EQ(m2(0, 1), 2.0);
	ASSERT_EQ(m2(0, 2), 3.0);
	ASSERT_EQ(m2(1, 0), 4.0);
	ASSERT_EQ(m2(1, 1), 5.0);
	ASSERT_EQ(m2(1, 2), 6.0);
	ASSERT_EQ(m2(2, 0), 7.0);
	ASSERT_EQ(m2(2, 1), 8.0);
	ASSERT_EQ(m2(2, 2), 9.0);

#if CONFIG_ENABLE_ASSERTIONS
	ASSERT_THROW((Matrix<2, 2>({1})), Exception);
	ASSERT_THROW((Matrix<2, 2>({1, 1, 1, 1, 1})), Exception);
#endif
}


TEST(Linal, MatrixAccess) {
	Matrix<3, 3> m({1.0, 2.0, 3.0,
	                4.0, 5.0, 6.0,
	                7.0, 8.0, 9.0});

	ASSERT_EQ(m(0, 0), 1.0);
	ASSERT_EQ(m(0, 1), 2.0);
	ASSERT_EQ(m(0, 2), 3.0);
	ASSERT_EQ(m(1, 0), 4.0);
	ASSERT_EQ(m(1, 1), 5.0);
	ASSERT_EQ(m(1, 2), 6.0);
	ASSERT_EQ(m(2, 0), 7.0);
	ASSERT_EQ(m(2, 1), 8.0);
	ASSERT_EQ(m(2, 2), 9.0);

	const auto& m_ref = m;

	ASSERT_EQ(m_ref(0, 2), 3.0);

	m(0, 2) = 10.0;
	ASSERT_EQ(m(0, 2), 10.0);


	SymmetricMatrix<3> sm({1.0, 2.0, 3.0,
	                            4.0, 5.0,
	                                 6.0});

	ASSERT_EQ(sm(0, 0), 1.0);
	ASSERT_EQ(sm(0, 1), 2.0);
	ASSERT_EQ(sm(0, 2), 3.0);
	ASSERT_EQ(sm(1, 0), 2.0);
	ASSERT_EQ(sm(1, 1), 4.0);
	ASSERT_EQ(sm(1, 2), 5.0);
	ASSERT_EQ(sm(2, 0), 3.0);
	ASSERT_EQ(sm(2, 1), 5.0);
	ASSERT_EQ(sm(2, 2), 6.0);
	
	const DiagonalMatrix<3> dm = {1, 2, 3};
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			ASSERT_EQ((i == j) * (i + 1), dm(i, j));
		}
	}
}


TEST(Linal, equalityInequality) {
	Matrix<2, 2> A = {1,  2, 3, 4};
	Matrix<2, 2> B = {1,  2, 3, 5};
	Matrix<2, 2> C = {1, -2, 9, 4};
	Matrix<2, 2> D = {0,  9, 9, 9};
		
	ASSERT_TRUE(A == A);
	ASSERT_TRUE(A != B);
	ASSERT_TRUE(A < B);
	ASSERT_TRUE(C < A);
	ASSERT_TRUE(C < B);
	ASSERT_TRUE(D < C);
	ASSERT_FALSE(B < A);
	ASSERT_FALSE(B < B);
	ASSERT_FALSE(A != A);
}


TEST(Linal, MatrixAdd) {
	Matrix<2, 2> m1({1.0, 2.0,
	                 3.0, 4.0});

	Matrix<2, 2> m2({4.0, 3.0,
	                 2.0, 1.0});
	
	SymmetricMatrix<2> s = {
			6, 7,
			   8,
	};

	auto m3 = m1 + m2;

	ASSERT_EQ(m3(0, 0), 5.0);
	ASSERT_EQ(m3(0, 1), 5.0);
	ASSERT_EQ(m3(1, 0), 5.0);
	ASSERT_EQ(m3(1, 1), 5.0);
	
	ASSERT_EQ(Matrix22({7, 9, 10, 12}), m1 + s);
	ASSERT_EQ(2 * s, s + s);
	
	MatrixInt<2, 2> mi = {1, 1, 2, 2};
	Matrix<2, 2> mr = {0.5, 1.5, 0, 2.2};
	
	auto res = mi + mr;
	ASSERT_EQ(res(0, 0), 1.5);
	ASSERT_EQ(res(0, 1), 2.5);
	ASSERT_EQ(res(1, 0), 2.0);
	ASSERT_EQ(res(1, 1), 4.2);
	
	auto res2 = mr + mi;
	ASSERT_EQ(res, res2);
}


TEST(Linal, MatrixSubtract) {
	Matrix<2, 2> m1({1.0, 2.0,
	                 3.0, 4.0});

	auto m2 = m1;


	auto m3 = m1 - m2;

	ASSERT_EQ(m3(0, 0), 0.0);
	ASSERT_EQ(m3(0, 1), 0.0);
	ASSERT_EQ(m3(1, 0), 0.0);
	ASSERT_EQ(m3(1, 1), 0.0);
}

TEST(Linal, MatrixNegative) {
	Matrix<2, 2> m1({1.0, 2.0,
	                 3.0, 4.0});

	auto m2 = -m1;

	ASSERT_EQ(m2(0, 0), -1.0);
	ASSERT_EQ(m2(0, 1), -2.0);
	ASSERT_EQ(m2(1, 0), -3.0);
	ASSERT_EQ(m2(1, 1), -4.0);
}


TEST(Linal, MatrixProduct) {
	Matrix<2, 2> m1({1.0, 2.0,
	                 3.0, 4.0});

	Matrix<2, 1> m2({5.0, 6.0});

	auto m3 = m1 * m2;

	ASSERT_EQ(m3(0, 0), 17.0);
	ASSERT_EQ(m3(1, 0), 39.0);
		
	MatrixInt<2, 3> i = {1, 0, 0,
	                     0, 0, -2};
	                     
	Matrix<3, 2> r = {1, 2,
	                  3, 4,
	                  5, 6};
	              
	SymmetricMatrix<2> s = {-2, 1, 0};
	
	DiagonalMatrix<2> d = {-5, 5};
	
	auto ir = i * r;
	ASSERT_EQ(Matrix22({1, 2, -10, -12}), ir);
	
	auto irs = ir * s;
	ASSERT_EQ(Matrix22({0, 1, 8, -10}), irs);
	
	auto irsd = irs * d;
	ASSERT_EQ(Matrix22({0, 5, -40, -50}), irsd);


	Matrix22 m4({1.0, 2.0,
	             3.0, 4.0});
	auto m5 = m4;
	auto m6 = m4 * m5;
	ASSERT_EQ(m6(0, 0), 7.0);
	ASSERT_EQ(m6(0, 1), 10.0);
	ASSERT_EQ(m6(1, 0), 15.0);
	ASSERT_EQ(m6(1, 1), 22.0);
}


TEST(Linal, DiagonalMatrix) {
	DiagonalMatrix<9> d = {1, 2, 3, 4, 5, 6, 7, 8, 9};
	ASSERT_EQ(Vector<9>({1, 2, 3, 4, 5, 6, 7, 8, 9}), diag(d));
	DiagonalMatrix<9> d2 = DiagonalMatrix<9>::Ones() * 2;
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 9; j++) {
			ASSERT_EQ((i == j) * 2 * (i + 1), (d*d2)(i, j));
			ASSERT_EQ((d*d2)(i, j), (d2*d)(i, j));
		}
	}
	SymmetricMatrix<3> s = {1, 2, 3, 4, 5, 6};
	ASSERT_EQ((DiagonalMatrix<3>({1, 4, 6})), Diag(s));
	Matrix33 m = {1, 2, 3, 4, 5, 6, 7, 8, 9};
	ASSERT_EQ((DiagonalMatrix<3>({1, 5, 9})), Diag(m));
}


TEST(Linal, MatrixEqual) {
	Matrix<2, 3> m1({1.0, 2.0, 3.0,
	                 4.0, 5.0, 6.0});

	auto m2 = m1;

	auto m3 = m1 * 2.0;

	ASSERT_EQ(m1, m2);
	ASSERT_NE(m1, m3);
	
	m1(0, 0) -= EQUALITY_TOLERANCE / 2;
	ASSERT_NE(m1, m2);
	ASSERT_TRUE(approximatelyEqual(m1, m2));
	m1(0, 0) -= EQUALITY_TOLERANCE;
	ASSERT_TRUE(approximatelyEqual(m1, m2, 2 * EQUALITY_TOLERANCE));
}


TEST(Linal, MatrixScalarMultiplication) {
	Matrix22 m({1.0, 2.0,
	            3.0, 4.0});

	m = 1.0 * m * 2.0;

	ASSERT_EQ(m(0, 0), 2.0);
	ASSERT_EQ(m(0, 1), 4.0);
	ASSERT_EQ(m(1, 0), 6.0);
	ASSERT_EQ(m(1, 1), 8.0);
}


TEST(Linal, MatrixScalarDivision) {
	Matrix22 m({2.0, 4.0,
	            6.0, 8.0});

	m = m / 2.0;

	ASSERT_EQ(m(0, 0), 1.0);
	ASSERT_EQ(m(0, 1), 2.0);
	ASSERT_EQ(m(1, 0), 3.0);
	ASSERT_EQ(m(1, 1), 4.0);
}


TEST(Linal, MatrixTranspose) {
	Matrix<2, 3> m1({1.0, 2.0, 3.0,
	                 4.0, 5.0, 6.0});

	auto m2 = transpose(m1);

	ASSERT_EQ(m2(0, 0), 1.0);
	ASSERT_EQ(m2(1, 0), 2.0);
	ASSERT_EQ(m2(2, 0), 3.0);
	ASSERT_EQ(m2(0, 1), 4.0);
	ASSERT_EQ(m2(1, 1), 5.0);
	ASSERT_EQ(m2(2, 1), 6.0);

	auto m3 = transpose(m2);

	ASSERT_EQ(m1, m3);
}


TEST(Linal, MatrixTransposeInplace) {
	Matrix<2, 2> m1({1.0, 2.0,
	                 3.0, 4.0});

	auto m2 = m1;

	m1.transposeInplace();

	ASSERT_EQ(m1(0, 0), 1.0);
	ASSERT_EQ(m1(1, 0), 2.0);
	ASSERT_EQ(m1(0, 1), 3.0);
	ASSERT_EQ(m1(1, 1), 4.0);

	m1.transposeInplace();

	ASSERT_EQ(m1, m2);
}


TEST(Linal, determinant) {
	Matrix33 m3({1.0, 2.0, 3.0,
	             3.0, 2.0, 1.0,
	             5.0, 7.0, 11.0});
	ASSERT_EQ(determinant(m3), -8.0);
	
	Matrix22 m2 = {1, 2, 3, 4};
	ASSERT_EQ(determinant(m2), -2);
}


TEST(Linal, Real3Construct) {
	Real3 m1({1.0, 2.0, 3.0});

	ASSERT_EQ(m1(0), 1.0);
	ASSERT_EQ(m1(1), 2.0);
	ASSERT_EQ(m1(2), 3.0);
	
	auto m2 = m1;

	ASSERT_EQ(m2(0), 1.0);
	ASSERT_EQ(m2(1), 2.0);
	ASSERT_EQ(m2(2), 3.0);
}


TEST(Linal, VectorLength) {
	Real3 v({0.0, 3.0, 4.0});

	ASSERT_EQ(length(v), 5.0);
}


TEST(Linal, VectorDotProduct) {
	SymmetricMatrix<3> H  = SymmetricMatrix<3>::Identity();
	SymmetricMatrix<3> H2 = {2, -1,  0,
	                             2, -1,
	                                 2};
	Real3 v1({1.0, 2.0, 3.0});
	Real3 v2({-2.0, 4.0, -2.0});
	ASSERT_EQ(dotProduct(v1, v2), 0.0);
	ASSERT_EQ(dotProduct(v1, H, v2), 0.0);
	ASSERT_EQ(dotProduct(v1, H2, v2), -8.0);
	
	Real3 v3({2, 2, 5});
	ASSERT_EQ(dotProduct(v1, v3), 21);
	ASSERT_EQ(dotProduct(v3, v1), 21);
	ASSERT_EQ(dotProduct(v3, H2, v1), 20);
	ASSERT_EQ(dotProduct(v1, H2, v3), 20);
}


TEST(Linal, VectorCrossProduct) {
	Real3 v1({1.0, 0.0, 0.0});

	Real3 v2({0.0, 1.0, 0.0});

	Real3 v3({0.0, 0.0, 1.0});

	Real3 z({0.0, 0.0, 0.0});

	ASSERT_EQ(crossProduct(v1, v2), v3);
	ASSERT_EQ(crossProduct(v2, v3), v1);
	ASSERT_EQ(crossProduct(v3, v1), v2);

	ASSERT_EQ(crossProduct(v2, v1), -v3);
	ASSERT_EQ(crossProduct(v3, v2), -v1);
	ASSERT_EQ(crossProduct(v1, v3), -v2);

	ASSERT_EQ(crossProduct(v1, v1), z);
	ASSERT_EQ(crossProduct(v2, v2), z);
	ASSERT_EQ(crossProduct(v3, v3), z);
	
	ASSERT_EQ( 1, crossProduct(Real2({1, 0}), Real2({1, 1})));
	ASSERT_EQ(-1, crossProduct(Real2({1, 1}), Real2({1, 0})));
	
	ASSERT_EQ(crossProduct(Real3({12, 5, 0}), Real3({-10, 1, 0}))(2),
	          crossProduct(Real2({12, 5}),    Real2({-10, 1})));
}


TEST(Linal, VectorNormalize) {
	Real3 v1({2.0, 0.0, 0.0});

	auto v = normalize(v1);

	ASSERT_NEAR(v(0), 1.0, 1e-5);
	ASSERT_EQ(v(1), 0.0);
	ASSERT_EQ(v(2), 0.0);

	Real3 v2({1.0, 2.0, 3.0});

	v = normalize(v2);

	ASSERT_NEAR(length(v), 1.0, 1e-5);

#if CONFIG_ENABLE_ASSERTIONS
	Real3 v3({0.0, 0.0, 0.0});

	ASSERT_THROW(normalize(v3), Exception);
#endif
}


TEST(Linal, RotationMatrix) {
	Real3 x({1.0, 0.0, 0.0});

	Real3 y({0.0, 1.0, 0.0});

	Real3 z({0.0, 0.0, 1.0});

	ASSERT_TRUE(approximatelyEqual(getXRotationMatrix(M_PI / 2) * y, -z)) 
		<< getXRotationMatrix(M_PI / 2) * y << "vs" << -z;
	ASSERT_TRUE(approximatelyEqual(getXRotationMatrix(M_PI / 2) * z, y));
	ASSERT_TRUE(approximatelyEqual(getXRotationMatrix(-M_PI / 2) * y, z));
	ASSERT_TRUE(approximatelyEqual(getXRotationMatrix(-M_PI / 2) * z, -y));

	ASSERT_TRUE(approximatelyEqual(getYRotationMatrix(M_PI / 2) * x, z));
	ASSERT_TRUE(approximatelyEqual(getYRotationMatrix(M_PI / 2) * z, -x));
	ASSERT_TRUE(approximatelyEqual(getYRotationMatrix(-M_PI / 2) * z, x));
	ASSERT_TRUE(approximatelyEqual(getYRotationMatrix(-M_PI / 2) * x, -z));

	ASSERT_TRUE(approximatelyEqual(getZRotationMatrix(M_PI / 2) * y, x));
	ASSERT_TRUE(approximatelyEqual(getZRotationMatrix(M_PI / 2) * x, -y));
	ASSERT_TRUE(approximatelyEqual(getZRotationMatrix(-M_PI / 2) * y, -x));
	ASSERT_TRUE(approximatelyEqual(getZRotationMatrix(-M_PI / 2) * x, y));
}


TEST(Linal, MatrixMatrixMultiplication) {
	Matrix<5, 5> A; Matrix<5, 5> B; Matrix<5, 5> C; // C = A * B

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

	Matrix<5, 5> AB = A * B;
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			ASSERT_NEAR(AB(i, j), C(i, j), EQUALITY_TOLERANCE);
		}
	}
}


TEST(Linal, MatrixVectorMultiplication) {
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


TEST(Linal, Trace) {
	Matrix<5, 5> A;
	A(0, 0) = 12; A(1, 1) = 56.333; A(2, 2) = 1; A(3, 3) = 0; A(4, 4) = -34.0022;
	ASSERT_NEAR(trace(A), 35.3308, EQUALITY_TOLERANCE);
}


TEST(Linal, setColumn) {
	Matrix<15, 33> matrix;
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


TEST(Linal, modifyAdd) {
	Matrix22 a = {
			1, 2, 
			3, 4,
	};
	
	SymmetricMatrix<2> s = {
			1, 2,
			   3,
	};
	
	a += s;
	ASSERT_EQ(Matrix22({2, 4, 5, 7}), a);
	
	a -= s;
	ASSERT_EQ(Matrix22({1, 2, 3, 4}), a);	
}


TEST(Linal, getColumn) {
	Matrix<22, 13> matrix;
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


TEST(Linal, setRow) {
	Matrix<15, 33> matrix;
	for (int i = 0; i < matrix.M; i++) {
		Vector<33> vector;
		for (int j = 0; j < matrix.N; j++) {
			vector(j) = i;
		}
		matrix.setRow(i, vector);
	}

	for (int i = 0; i < matrix.M; i++) {
		for (int j = 0; j < matrix.N; j++) {
			ASSERT_EQ(matrix(i, j), i);
		}
	}
}


TEST(Linal, getRow) {
	Matrix<22, 13> matrix;
	for (int i = 0; i < matrix.M; i++) {
		for (int j = 0; j < matrix.N; j++) {
			matrix(i, j) = i;
		}
	}

	for (int k = 0; k < matrix.M; k++) {
		Vector<13> string = matrix.getRow(k);
		for (int i = 0; i < string.M; i++) {
			ASSERT_EQ(string(i), k);
		}
	}
}


TEST(Linal, diagonalMultiply) {
	Matrix<9, 9> A, B;
	Utils::seedRand();

	for (int l = 0; l < LINAL_TEST_NUMBER_ITERATIONS; l++) {

		for (int i = 0; i < A.M; i++) {
			for (int j = 0; j < A.N; j++) {
				A(i, j) = Utils::randomReal(-1e+6, 1e+6);
				B(i, j) = Utils::randomReal(-1e+6, 1e+6);
			}
		}

		Matrix<9, 9> C = A * B;
		Vector<9> d = diagonalMultiply(A, B);

		for (int k = 0; k < d.M; k++) {
			ASSERT_EQ(C(k, k), d(k));
		}

	}
}


TEST(Linal, diag) {
	Matrix33 m = {1, 2, 3, 4, 5, 6, 7, 8, 9};
	ASSERT_EQ(Real3({1, 5, 9}), diag(m));
}


TEST(Linal, VectorOperators) {
	Vector<7> vector;
	for (int i = 0; i < vector.M; i++) {
		vector(i) = (real) i - 3;
	}

	vector = vector * 2.0;
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
	
	vector -= vector;
	ASSERT_EQ(Vector<7>::Zeros(), vector);
}


TEST(Linal, VectorInt) {
	VectorInt<3> i = {1, 3, 4};
	Vector<3> r = {0.1, 0.3, 0.7};

	VectorInt<3> plainIntMultpl = plainMultiply(i, r);
	ASSERT_EQ(0, plainIntMultpl(0));
	ASSERT_EQ(0, plainIntMultpl(1));
	ASSERT_EQ(2, plainIntMultpl(2));

	Vector<3> plainRealMultpl = plainMultiply(i, r);
	// attention: NEAR!!!
	ASSERT_NEAR(0.1, plainRealMultpl(0), EQUALITY_TOLERANCE);
	ASSERT_NEAR(0.9, plainRealMultpl(1), EQUALITY_TOLERANCE);
	ASSERT_NEAR(2.8, plainRealMultpl(2), EQUALITY_TOLERANCE);
}


TEST(Linal, plainDivision) {
	VectorInt<3> i = {1, 3, 4};
	Vector<3> r = {0.1, 0.3, 0.7};

	Vector<3> q = plainDivision(r, i);
	ASSERT_NEAR(0.1, q(0), EQUALITY_TOLERANCE);
	ASSERT_NEAR(0.1, q(1), EQUALITY_TOLERANCE);
	ASSERT_NEAR(0.7 / 4, q(2), EQUALITY_TOLERANCE);

	Vector<1> a = {1};
	Vector<1> b = {0};
	ASSERT_NEAR(std::numeric_limits<real>::max(),
	            plainDivision(a, b)(0), EQUALITY_TOLERANCE);
	ASSERT_NEAR(-std::numeric_limits<real>::max(),
	            plainDivision(-a, b) (0), EQUALITY_TOLERANCE);
}


TEST(Linal, directProduct) {
	VectorInt<4> i = {1, 2, 3, 4};
	ASSERT_EQ(24, directProduct(i));
	
	Vector<8> r = {1, 2, 3, 4, 5, 6, 7, 0.5};
	ASSERT_EQ(2520, directProduct(r));
}


TEST(Linal, perpendicularClockwise) {
	ASSERT_EQ(Real2({0, -1}), perpendicularClockwise(Real2({1, 0})));
	ASSERT_EQ(Real2({5, 0}), perpendicularClockwise(Real2({0, 5})));
	ASSERT_EQ(Real2({3, -2}), perpendicularClockwise(Real2({2, 3})));
}


TEST(Linal, solveLinearSystem) {
	Matrix<1, 1> A1 = {5};
	Vector<1> b1 = {5};
	auto x1 = solveLinearSystem(A1, b1);
	ASSERT_EQ(Real1({1}), x1);

	Matrix<2, 2> A2 = {1, 2,
	                   3, 4};
	Vector<2> b2 = {5, 6};
	auto x2 = solveLinearSystem(A2, b2);
	ASSERT_EQ(Real2({-4, 4.5}), x2);
	
	Matrix<3, 3> A3 = {1, 3, -2,
	                   3, 5, 6,
	                   2, 4, 3};
	Vector<3> b3 = {5, 7, 8};
	auto x3 = solveLinearSystem(A3, b3);
	ASSERT_EQ(Real3({-15, 8, 2}), x3);
}


TEST(Linal, createLocalBasis) {
	auto b1 = createLocalBasis(Real1({1}));
	ASSERT_EQ(Matrix11({1}), b1);
	auto b1T = createLocalBasisTranspose(Real1({1}));
	ASSERT_EQ(Matrix11({1}), b1T);
	
	
	auto b2 = createLocalBasis(Real2({0, 1}));
	ASSERT_EQ(Matrix22({1, 0,
	                    0, 1}), b2);
	
	b2 = createLocalBasis(Real2({1, 0}));
	ASSERT_EQ(Matrix22({0, 1,
	                   -1, 0}), b2);
	
	auto b2T = createLocalBasisTranspose(Real2({1, 0}));
	ASSERT_EQ(transpose(b2), b2T);
	
	
	auto b3 = createLocalBasis(Real3({0, 0, 1}));
	ASSERT_EQ(Matrix33({1, 0, 0,
	                    0, 1, 0,
	                    0, 0, 1}), b3);

	b3 = createLocalBasis(Real3({0, 1, 0}));
	ASSERT_EQ(Matrix33({1,  0, 0,
	                    0,  0, 1,
	                    0, -1, 0}), b3);	

	b3 = createLocalBasis(Real3({1, 0, 0}));
	ASSERT_EQ(Matrix33({0,  0, 1,
	                   -1,  0, 0,
	                    0, -1, 0}), b3);
	
	auto b3T = createLocalBasisTranspose(Real3({1, 0, 0}));
	ASSERT_EQ(transpose(b3), b3T);
}


TEST(Linal, cross_verification) {
	Matrix<3, 3> A, B, C, D1;
	DiagonalMatrix<3> D;
	Vector<3> f;
	Utils::seedRand();

	for (int l = 0; l < LINAL_TEST_NUMBER_ITERATIONS; l++) {

		for (int i = 0; i < A.M; i++) {
			for (int j = 0; j < A.N; j++) {
				A(i, j) = Utils::randomReal(-1e+6, 1e+6);
				B(i, j) = Utils::randomReal(-1e+6, 1e+6);
				D1(i, j) = 0;
			}
			D1(i, i) = D(i) = Utils::randomReal(-1e+6, 1e+6);
			f(i) = Utils::randomReal(-1e+6, 1e+6);
		}
		
		// plus - minus
		C = A + B;
		ASSERT_TRUE(approximatelyEqual(A, C - B));
		ASSERT_TRUE(approximatelyEqual(B, C - A));
		C -= B;
		ASSERT_TRUE(approximatelyEqual(A, C));
		C -= A;
		ASSERT_TRUE(approximatelyEqual(zeros(C), C, 1000 * EQUALITY_TOLERANCE));
		C = 2.0 * A;
		C = C - 3.0 * A;
		ASSERT_TRUE(approximatelyEqual(- A, C));
		
		
		// multiplication
		C = A * B;
		ASSERT_EQ(transpose(B) * transpose(A), transpose(C));
		ASSERT_EQ(transpose(A) * B, transposeMultiply(A, B));
		ASSERT_NEAR(determinant(B) * determinant(A), 
		            determinant(C), EQUALITY_TOLERANCE  * fabs(determinant(C)));
		auto C1 = GslUtils::invert(C);
		ASSERT_TRUE(approximatelyEqual(GslUtils::invert(B) * GslUtils::invert(A), C1));
		ASSERT_TRUE(approximatelyEqual(identity(C), C * C1, 1000 * EQUALITY_TOLERANCE));
		ASSERT_TRUE(approximatelyEqual(identity(C), C1 * C, 1000 * EQUALITY_TOLERANCE));
		ASSERT_EQ(A * D1, A * D);
		ASSERT_EQ(D1 * A, D * A);
		
		// SLE
		auto x = solveLinearSystem(A, f);
		ASSERT_TRUE(approximatelyEqual(f, A * x)) << "Ax - f = " << A * x - f;
		ASSERT_TRUE(approximatelyEqual(f, solveLinearSystem(A, A * f)));
		// SLE + least squares
		ASSERT_TRUE(approximatelyEqual(solveLinearSystem(A, f), 
		                               linearLeastSquares(A, f),
		                               1000 * EQUALITY_TOLERANCE))
			<< solveLinearSystem(A, f) << "vs" << linearLeastSquares(A, f);
	}
}


TEST(Linal, SLE_with_vectors_at_right_side) {
	Matrix22 A = {1, 0,
	              0, 1};
	Vector<3> b1 = {1, 2, 3};
	Vector<3> b2 = {4, 5, 6};
	
	VECTOR<2, Vector<3>> b = {b1, b2};
	
	auto x = solveLinearSystem(A, b);
	
	ASSERT_EQ(b1, x(0));
	ASSERT_EQ(b2, x(1));
	
	A = {1, 2,
	     1, 0};
	b = {{1, 0, 1}, {0, 0, 2}};
	x = solveLinearSystem(A, b);
	ASSERT_EQ(Vector<3>({0,   0,  2  }), x(0));
	ASSERT_EQ(Vector<3>({0.5, 0, -0.5}), x(1));
}


TEST(Linal, linearLeastSquares) {
	real a1 = 1, a2 = 2, b1 = 3, b2 = 4;
	Matrix<2, 1> A = {a1, a2};
	Vector<2> b = {b1, b2};
	
	auto x = linearLeastSquares(A, b);
	ASSERT_EQ((a1*b1 + a2*b2) / (a1*a1 + a2*a2), x(0));
	
	Matrix<5, 2> A5 = {1, 1,
	                   2, 0,
	                   3, 0,
	                   4, 0,
	                   5, 0};
	Vector<5> b5 = {1, 2, 3, 4, 5};
	ASSERT_EQ(Vector<2>({1, 0}), linearLeastSquares(A5, b5));
	ASSERT_EQ(Vector<2>({1, 0}), linearLeastSquares(A5, b5, DiagonalMatrix<5>::Identity()));
	
	A = {1, 2};
	b = {1, 1};
	DiagonalMatrix<2> W = {7, 0};
	ASSERT_EQ(1, linearLeastSquares(A, b, W)(0));
	W = {0, 1};
	ASSERT_EQ(0.5, linearLeastSquares(A, b, W)(0));
	W = {1, 2};
	ASSERT_EQ(5.0 / 9.0, linearLeastSquares(A, b, W)(0));
}


TEST(Linal, linesIntersection2D) {
	ASSERT_EQ(Real2({1, 1}), 
			linesIntersection(Real2({1, 0}),  Real2({1, 2}), Real2({0, 1}), Real2({2, 1})));
	ASSERT_EQ(Real2({0, 0}), 
			linesIntersection(Real2({-1, -1}), Real2({1, 1}), Real2({-1, 1}), Real2({1, -1})));
	ASSERT_EQ(Real2({1, 1}), 
			linesIntersection(Real2({-1, -1}), Real2({-2, -2}), Real2({1, 0}), Real2({1, -3})));
}


TEST(Linal, linesIntersection3D) {
	ASSERT_EQ(Real3({1, 1, 0}), 
			linesIntersection(Real3({1, 0, 0}),  Real3({1, 2, 0}),
					Real3({0, 1, 0}), Real3({2, 1, 0})));
	ASSERT_EQ(Real3({0, 0, 5}), 
			linesIntersection(Real3({-1, -1, 5}), Real3({1, 1, 5}),
					Real3({-1, 1, 5}), Real3({1, -1, 5})));
	ASSERT_EQ(Real3({0, 0, 0}), 
			linesIntersection(Real3({-1, -1, -1}), Real3({1, 1, 1}),
					Real3({-1, 0, 0}), Real3({1, 0, 0})));
}


TEST(Linal, lineWithFlatIntersection) {
	ASSERT_EQ(Real3({0, 0, 0}), lineWithFlatIntersection(
			{1, 0, 0}, {0, 1, 0}, {0, 0, 0}, {1, 1, 1}, {2, 2, 2}));
	ASSERT_EQ(Real3({1, 1, 1}), lineWithFlatIntersection(
			{3, 0, 0}, {0, 3, 0}, {0, 0, 3}, {-1, -1, -1}, {5, 5, 5}));
	ASSERT_EQ(Real3({8, 2, 3}), lineWithFlatIntersection(
			{8, 2, 3}, {4, 5, 6}, {9, 8, 7}, {8, 2, 3}, {2, 2, 2}));
}


TEST(Linal, barycentric1D) {
	ASSERT_EQ(Real2({ 1,  0}), barycentricCoordinates(Real1({3}), Real1({4}), Real1({3})));
	ASSERT_EQ(Real2({ 0,  1}), barycentricCoordinates(Real1({3}), Real1({4}), Real1({4})));
	ASSERT_EQ(Real2({-1,  2}), barycentricCoordinates(Real1({3}), Real1({4}), Real1({5})));
	ASSERT_EQ(Real2({ 2, -1}), barycentricCoordinates(Real1({3}), Real1({4}), Real1({2})));
	ASSERT_TRUE(linal::approximatelyEqual(Real2({ 0.3, 0.7 }), 
			barycentricCoordinates(Real1({3}), Real1({4}), Real1({3.7}))));
}


TEST(Linal, barycentricSegmentIn2D) {
	Real2 a = {0, 0}, b = {0, 1}, c = {4, 7};
	ASSERT_EQ(Real2({ 1,  0}), barycentricCoordinates(a, b, a));
	ASSERT_EQ(Real2({ 0,  1}), barycentricCoordinates(a, b, b));
	ASSERT_EQ(Real2({ 0.5, 0.5}), barycentricCoordinates(a, b, (a + b) / 2));
	ASSERT_EQ(Real2({ 0.5, 0.5}), barycentricCoordinates(c, b, (c + b) / 2));
	ASSERT_EQ(Real2({ -2, 3}), barycentricCoordinates(a, c, 3 * c));
}


TEST(Linal, barycentric2D) {
	Real2 a = {0, 0};
	Real2 b = {1, 0};
	Real2 c = {0, 1};
	
	ASSERT_EQ(Real3({1, 0, 0}), barycentricCoordinates(a, b, c, a));
	ASSERT_EQ(Real3({0, 1, 0}), barycentricCoordinates(a, b, c, b));
	ASSERT_EQ(Real3({0, 0, 1}), barycentricCoordinates(a, b, c, c));

	ASSERT_EQ(Real3({0, 0.5, 0.5}), barycentricCoordinates(a, b, c, (b+c)/2.0));
	ASSERT_TRUE(approximatelyEqual(Real3({1.0 / 3, 1.0 / 3, 1.0 / 3}),
			barycentricCoordinates(a, b, c, (a+b+c)/3.0)));

	// affine invariance	
	Matrix22 S = {3, 5,
	             -2, 10};
	Real2 shift = {5, -6};
	
	a = S * a + shift;
	b = S * b + shift;
	c = S * c + shift;
	
	ASSERT_EQ(Real3({1, 0, 0}), barycentricCoordinates(a, b, c, a));
	ASSERT_EQ(Real3({0, 1, 0}), barycentricCoordinates(a, b, c, b));
	ASSERT_EQ(Real3({0, 0, 1}), barycentricCoordinates(a, b, c, c));

	ASSERT_EQ(Real3({0, 0.5, 0.5}), barycentricCoordinates(a, b, c, (b+c)/2.0));
	ASSERT_TRUE(approximatelyEqual(Real3({1.0 / 3, 1.0 / 3, 1.0 / 3}),
			barycentricCoordinates(a, b, c, (a+b+c)/3.0)));
}


TEST(Linal, barycentricTriangleIn3D) {
	Real3 a = {0, 0, 1};
	Real3 b = {1, 0, 0};
	Real3 c = {0, 1, 0};
	
	ASSERT_TRUE(approximatelyEqual(Real3({0, 1, 0}),
			barycentricCoordinates(a, b, c, b)));
	ASSERT_TRUE(approximatelyEqual(Real3({0, 0, 1}),
			barycentricCoordinates(a, b, c, c)));
	ASSERT_TRUE(approximatelyEqual(Real3({1, 0, 0}),
			barycentricCoordinates(a, b, c, a)));
	ASSERT_TRUE(approximatelyEqual(Real3({0, 0.5, 0.5}),
			barycentricCoordinates(a, b, c, (b + c) / 2.0)));
	ASSERT_TRUE(approximatelyEqual(Real3({1.0 / 3, 1.0 / 3, 1.0 / 3}),
			barycentricCoordinates(a, b, c, (a + b + c) / 3.0)));
	
	a = {-1.587519073486328, 0.00633697509765625, 143.9256103515625};
	b = {-1.584312438964844, 0.0194580078125, 143.9030517578125};
	c = {-1.575632476806641, 0.0271484375, 143.9210571289063};
	Real3 q = {-1.581665132409098, 0.0211456298828125, 143.9104125976563};
	ASSERT_TRUE(isDegenerate(a, b, c, q, 7e-11));
	Real3 lambdaExpt = {0.063996977831142207, 0.60737008320664498, 0.32863293896221279};
	Real3 lambdaCalc = barycentricCoordinates(a, b, c, q);
	ASSERT_LT(linal::length(lambdaCalc - lambdaExpt), EQUALITY_TOLERANCE)
			<< "lambdaCalc: " << lambdaCalc << "lambdaExpt: " << lambdaExpt;
	Real3 qCalc = lambdaCalc(0) * a + lambdaCalc(1) * b + lambdaCalc(2) * c;
	ASSERT_LT(linal::length(qCalc - q), EQUALITY_TOLERANCE)
			<< "qCalc: " << qCalc << "q: " << q;
}


TEST(Linal, barycentric3D) {
	Real3 a = {0, 0, 0};
	Real3 b = {1, 0, 0};
	Real3 c = {0, 1, 0};
	Real3 d = {0, 0, 1};
	
	ASSERT_EQ(Real4({1, 0, 0, 0}), barycentricCoordinates(a, b, c, d, a));
	ASSERT_EQ(Real4({0, 1, 0, 0}), barycentricCoordinates(a, b, c, d, b));
	ASSERT_EQ(Real4({0, 0, 1, 0}), barycentricCoordinates(a, b, c, d, c));
	ASSERT_EQ(Real4({0, 0, 0, 1}), barycentricCoordinates(a, b, c, d, d));

	ASSERT_EQ(Real4({0, 0.5, 0.5, 0}), barycentricCoordinates(a, b, c, d, (b+c)/2.0));
	ASSERT_EQ(Real4({0.25, 0.25, 0.25, 0.25}), barycentricCoordinates(a, b, c, d, (a+b+c+d)/4));

	// affine invariance
	Matrix33 S = {3, 5, 2,
	             -2, 10, 0,
	              1, 0, -6};
	Real3 shift = {5, -6, 1};
	
	a = S * a + shift;
	b = S * b + shift;
	c = S * c + shift;
	d = S * d + shift;
	
	ASSERT_EQ(Real4({1, 0, 0, 0}), barycentricCoordinates(a, b, c, d, a));
	ASSERT_EQ(Real4({0, 1, 0, 0}), barycentricCoordinates(a, b, c, d, b));
	ASSERT_EQ(Real4({0, 0, 1, 0}), barycentricCoordinates(a, b, c, d, c));
	ASSERT_EQ(Real4({0, 0, 0, 1}), barycentricCoordinates(a, b, c, d, d));

	ASSERT_EQ(Real4({0, 0.5, 0.5, 0}), barycentricCoordinates(a, b, c, d, (b+c)/2.0));
	ASSERT_EQ(Real4({0.25, 0.25, 0.25, 0.25}), barycentricCoordinates(a, b, c, d, (a+b+c+d)/4));
}


TEST(Linal, isPerpendicular) {
	ASSERT_TRUE (isPerpendicular(Real2({ 0, 1}), Real2({1, 0})));
	ASSERT_FALSE(isPerpendicular(Real2({ 0, 1}), Real2({1, 1})));
	ASSERT_TRUE (isPerpendicular(Real3({0, 1, 0}), Real3({1, 0, 1})));
	ASSERT_FALSE(isPerpendicular(Real3({1, 0, 1}), Real3({1, 0, 0})));
}


TEST(Linal, isDegenerate) {
	ASSERT_TRUE (isDegenerate(Real2({0, 0}), Real2({0, 0}), Real2({0, 0}), 0));
	ASSERT_TRUE (isDegenerate(Real2({0, 0}), Real2({0, 1}), Real2({0, 2}), 0));
	ASSERT_TRUE (isDegenerate(Real2({0, 0}), Real2({1, 1}), Real2({2, 2}), 0));
	ASSERT_TRUE (isDegenerate(Real2({0, 0}), Real2({0, 0}), Real2({1, 1}), 0));
	ASSERT_FALSE(isDegenerate(Real2({0, 0}), Real2({1, 0}), Real2({0, 1}), 0));
	ASSERT_FALSE(isDegenerate(Real2({0, 0}), Real2({-1, 1e-7}), Real2({1, 1e-7}), 0));
	ASSERT_TRUE (isDegenerate(Real2({0, 0}), Real2({-1, 1e-7}), Real2({1, 1e-7}), 2e-7));
	
	ASSERT_TRUE (isDegenerate(Real3({0, 0, 0}), Real3({0, 0, 1}), Real3({0, 0, 2}), 0));
	ASSERT_TRUE (isDegenerate(Real3({0, 0, 0}), Real3({1, 1, 1}), Real3({2, 2, 2}), 0));
	ASSERT_FALSE(isDegenerate(Real3({0, 5, 0}), Real3({0, 0, 1}), Real3({0, 0, 2}), 0));
	ASSERT_FALSE(isDegenerate(Real3({0, 0, 0}), Real3({1, 1, 1}), Real3({2, 2, 2 + 1e-7}), 0));
	ASSERT_TRUE (isDegenerate(Real3({0, 0, 0}), Real3({1, 1, 1}), Real3({2, 2, 2 + 1e-7}), 2e-7));
	
	ASSERT_TRUE (isDegenerate(Real3({0, 0, 0}), Real3({0, 0, 1}), Real3({0, 0, 2}), Real3({0, 2, 3}), 0));
	ASSERT_TRUE (isDegenerate(Real3({3, 0, 0}), Real3({0, 3, 0}), Real3({0, 0, 3}), Real3({1, 1, 1}), 0));
	ASSERT_FALSE(isDegenerate(Real3({3, 0, 0}), Real3({0, 3, 0}), Real3({0, 0, 3}), Real3({1, 1, 5}), 0));
	ASSERT_FALSE(isDegenerate(Real3({1e-7, 0, 0}), Real3({0, 0, 1}), Real3({0, 0, 2}), Real3({0, 2, 3}), 0));
	ASSERT_TRUE (isDegenerate(Real3({1e-7, 0, 0}), Real3({0, 0, 1}), Real3({0, 0, 2}), Real3({0, 2, 3}), 2e-7));
}


TEST(Linal, orientedArea) {
	ASSERT_EQ( 0, orientedArea({0, 0}, {0, 1}, {0, 2}));
	ASSERT_EQ( 0, orientedArea({0, 0}, {1, 1}, {2, 2}));
	ASSERT_EQ( 0, area(Real2({0, 0}), Real2({1, 1}), Real2({2, 2})));
	ASSERT_EQ(-2, orientedArea({0, 0}, {0, 2}, {2, 0}));
	ASSERT_EQ( 2, area(Real2({0, 0}), Real2({0, 2}), Real2({2, 0})));
	ASSERT_EQ( 2, orientedArea({0, 0}, {2, 0}, {0, 2}));
	ASSERT_EQ( 2, area(Real2({0, 0}), Real2({2, 0}), Real2({0, 2})));
	
	ASSERT_EQ( 2, area(Real3({0, 0, 0}), Real3({0, 2, 0}), Real3({2, 0, 0})));
	ASSERT_NEAR(sqrt(3) / 2, area(Real3({1, 0, 0}), Real3({0, 1, 0}), Real3({0, 0, 1})), EQUALITY_TOLERANCE);
}


TEST(Linal, orientedVolume) {
	ASSERT_EQ( 0, orientedVolume({0, 0, 0}, {0, 1, 0}, {0, 2, 0}, {1, 2, 0}));
	ASSERT_EQ( 0, volume({0, 0, 0}, {0, 1, 0}, {0, 2, 0}, {1, 2, 0}));
	ASSERT_EQ( 1.0 / 6, orientedVolume({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}));
	ASSERT_EQ( 1.0 / 6, volume({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}));
	ASSERT_EQ(-1.0 / 6, orientedVolume({0, 0, 0}, {0, 1, 0}, {1, 0, 0}, {0, 0, 1}));
	ASSERT_EQ( 1.0 / 6, volume({0, 0, 0}, {0, 1, 0}, {1, 0, 0}, {0, 0, 1}));
	ASSERT_EQ(-1.0 / 6, orientedVolume({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0,-1}));
	ASSERT_EQ( 1.0 / 6, volume({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0,-1}));
}


TEST(Linal, minimalHeight) {
	ASSERT_NEAR(1.0 / sqrt(2), minimalHeight({0, 0}, {1, 0}, {0, 1}), EQUALITY_TOLERANCE);
	ASSERT_NEAR(2, minimalHeight({0, 0}, {5, 2}, {8, 0}), EQUALITY_TOLERANCE);
	
	ASSERT_NEAR(1.0 / sqrt(3), minimalHeight({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}), EQUALITY_TOLERANCE);
	ASSERT_NEAR(2, minimalHeight({0, 0, 0}, {5, 0, 0}, {0, 8, 0}, {2, 2, 2}), EQUALITY_TOLERANCE);
}


TEST(Linal, oppositeFaceNormal) {
	ASSERT_EQ(Real3({1, 0, 0}), oppositeFaceNormal(
			{-1, 0, 0}, {0, 0, 0}, {0, 1, 0}, {0, 0, 1}));
	ASSERT_EQ(Real3({-1, 0, 0}), oppositeFaceNormal(
			{1, 0, 0}, {0, 0, 0}, {0, 1, 0}, {0, 0, 1}));
	ASSERT_EQ(normalize(Real3({1, 1, 1})), oppositeFaceNormal(
			{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}));
}


TEST(Linal, solveDegenerateLinearSystem) {
	ASSERT_EQ(Real3({0, 0, 1}), 
			solveDegenerateLinearSystem(Matrix33({1, 2, 0,
	                                              0, 0, 0,
	                                              1, 3, 0}), 1).at(0));
	ASSERT_EQ(Real3({0, 1, 0}), 
			solveDegenerateLinearSystem(Matrix33({1, 0, 0,
	                                              0, 0, 0,
	                                              0, 0, 0}), 2).at(0));
	ASSERT_EQ(Real3({0, 0, 1}), 
			solveDegenerateLinearSystem(Matrix33({1, 0, 0,
	                                              0, 0, 0,
	                                              0, 0, 0}), 2).at(1));
	ASSERT_EQ(Real3({1, -13.0 / 14, 4.0 / 14}), 
			solveDegenerateLinearSystem(Matrix33({ 1,    2,    3,
	                                               5,    6,    2,
	                                               1,    2,    3}), 1).at(0));
	ASSERT_EQ(Real3({-0.5, 1, 0}), 
			solveDegenerateLinearSystem(Matrix33({ 2,     1,    1,
	                                               2,     1,    1,
	                                               2,     1,    1,}), 2).at(0));
	ASSERT_EQ(Real3({-0.5, 0, 1}), 
			solveDegenerateLinearSystem(Matrix33({ 2,     1,    1,
	                                               2,     1,    1,
	                                               2,     1,    1,}), 2).at(1));
}


TEST(Linal, symmDirectProduct) {
	Real3 n0 = {1, 2, 3};
	Real3 n1 = {10, 100, 1000};
	Real2 n2 = {10, -10};
	
	ASSERT_EQ(Matrix33({10, 100, 1000,
	                    20, 200, 2000,
	                    30, 300, 3000,}), directProduct(n0, n1));
	
	ASSERT_EQ((Matrix<3, 2>({10, -10, 
	                         20, -20,
	                         30, -30,})), directProduct(n0, n2));
	
	ASSERT_EQ((Matrix<2, 3>({ 10,  20,  30,
	                         -10, -20, -30,})), directProduct(n2, n0));

	ASSERT_EQ(Matrix33({10,  60,   515,
	                    60,  200,  1150,
	                    515, 1150, 3000 }), symmDirectProduct(n0, n1));
	
	ASSERT_EQ(SymmetricMatrix<3>({10,  60,   515,
	                                   200,  1150,
	                                         3000 }), symmDirectProduct(n0, n1));
	
	ASSERT_EQ(directProduct(n0, n0), symmDirectProduct(n0, n0));
}



