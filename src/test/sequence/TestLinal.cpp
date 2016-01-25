#include <lib/linal/Linal.hpp>
#include <lib/linal/special/Matrix33.hpp>
#include <lib/linal/special/Matrix22.hpp>
#include <lib/linal/special/Vector3.hpp>
#include <lib/linal/special/VectorInt.hpp>
#include <lib/linal/special/RotationMatrix.hpp>

#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

const int LINAL_TEST_NUMBER_ITERATIONS = 1000;

using namespace gcm;
using namespace gcm::linal;

template<int M, int N>
class MatrixWrapper: public gcm::linal::Matrix<M, N>
{
    public:
        int _getIndex(int i, int j) const
        {
            return this->getIndex(i, j);
        }
};

TEST(Linal, ContainerIndex)
{
    MatrixWrapper<2, 8> m1;
    MatrixWrapper<8, 2> m2;

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 8; j++)
        {
            ASSERT_LT(m1._getIndex(i, j), 16);
            ASSERT_LT(m2._getIndex(j, i), 16);
        }
}

TEST(Linal, MatrixConstruct)
{
    gcm::linal::Matrix<3, 3> m1({
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0
    });

    ASSERT_EQ(m1.values[0], 1.0);
    ASSERT_EQ(m1.values[1], 2.0);
    ASSERT_EQ(m1.values[2], 3.0);
    ASSERT_EQ(m1.values[3], 4.0);
    ASSERT_EQ(m1.values[4], 5.0);
    ASSERT_EQ(m1.values[5], 6.0);
    ASSERT_EQ(m1.values[6], 7.0);
    ASSERT_EQ(m1.values[7], 8.0);
    ASSERT_EQ(m1.values[8], 9.0);


    gcm::linal::Matrix<3, 3, Matrix33Container> m2({
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0
    });
    
    ASSERT_EQ(m2.xx, 1.0);
    ASSERT_EQ(m2.xy, 2.0);
    ASSERT_EQ(m2.xz, 3.0);
    ASSERT_EQ(m2.yx, 4.0);
    ASSERT_EQ(m2.yy, 5.0);
    ASSERT_EQ(m2.yz, 6.0);
    ASSERT_EQ(m2.zx, 7.0);
    ASSERT_EQ(m2.zy, 8.0);
    ASSERT_EQ(m2.zz, 9.0);
}

TEST(Linal, MatrixAccess)
{
    gcm::linal::Matrix<3, 3> m({
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0
    });

    ASSERT_EQ(m(0, 0), 1.0);

    const auto& m_ref = m;

    ASSERT_EQ(m_ref(0, 0), 1.0);

    m(0, 0) = 10.0;
    ASSERT_EQ(m(0, 0), 10.0);
}


TEST(Linal, MatrixAdd)
{
    gcm::linal::Matrix<2, 2> m1({
        1.0, 2.0,
        3.0, 4.0
    });
    
    gcm::linal::Matrix<2, 2> m2({
        4.0, 3.0,
        2.0, 1.0
    });

    auto m3 = m1 + m2;

    ASSERT_EQ(m3(0, 0), 5.0);
    ASSERT_EQ(m3(0, 1), 5.0);
    ASSERT_EQ(m3(1, 0), 5.0);
    ASSERT_EQ(m3(1, 1), 5.0);
}


TEST(Linal, MatrixAddCustomContainer)
{
    gcm::linal::Matrix<2, 2, gcm::linal::Matrix22Container> m1({
        1.0, 2.0,
        3.0, 4.0
    });
    
    gcm::linal::Matrix<2, 2, Matrix22Container> m2({
        4.0, 3.0,
        2.0, 1.0
    });

    auto m3 = m1 + m2;

    ASSERT_EQ(m3.xx, 5.0);
    ASSERT_EQ(m3.xy, 5.0);
    ASSERT_EQ(m3.yx, 5.0);
    ASSERT_EQ(m3.yy, 5.0);
}

TEST(Linal, MatrixAssign)
{
    gcm::linal::Matrix<2, 2> m1({
        1.0, 2.0,
        3.0, 4.0
    });

    auto m2 = m1;

    ASSERT_EQ(m2(0, 0), 1.0);
    ASSERT_EQ(m2(0, 1), 2.0);
    ASSERT_EQ(m2(1, 0), 3.0);
    ASSERT_EQ(m2(1, 1), 4.0);
}

TEST(Linal, MatrixSubtract)
{
    gcm::linal::Matrix<2, 2> m1({
        1.0, 2.0,
        3.0, 4.0
    });

    auto m2 = m1;
    

    auto m3 = m1 - m2;

    ASSERT_EQ(m3(0, 0), 0.0);
    ASSERT_EQ(m3(0, 1), 0.0);
    ASSERT_EQ(m3(1, 0), 0.0);
    ASSERT_EQ(m3(1, 1), 0.0);
}

TEST(Linal, MatrixNegative)
{
    gcm::linal::Matrix<2, 2> m1({
        1.0, 2.0,
        3.0, 4.0
    });

    auto m2 = -m1;
    
    ASSERT_EQ(m2(0, 0), -1.0);
    ASSERT_EQ(m2(0, 1), -2.0);
    ASSERT_EQ(m2(1, 0), -3.0);
    ASSERT_EQ(m2(1, 1), -4.0);
}

TEST(Linal, MatrixSubtractCustomContainer)
{
    gcm::linal::Matrix<2, 2, Matrix22Container> m1({
        1.0, 2.0,
        3.0, 4.0
    });

    auto m2 = m1;

    auto m3 = m1 - m2;

    ASSERT_EQ(m3.xx, 0.0);
    ASSERT_EQ(m3.xy, 0.0);
    ASSERT_EQ(m3.yx, 0.0);
    ASSERT_EQ(m3.yy, 0.0);
}

TEST(Linal, MatrixProduct)
{
    gcm::linal::Matrix<2, 2> m1({
        1.0, 2.0,
        3.0, 4.0
    });

    gcm::linal::Matrix<2, 1> m2({
        5.0, 6.0
    });

    auto m3 = m1 * m2;

    ASSERT_EQ(m3(0, 0), 17.0);
    ASSERT_EQ(m3(1, 0), 39.0);
}

TEST(Linal, MatrixProductSquare)
{
    gcm::linal::Matrix<2, 2, Matrix22Container> m1({
        1.0, 2.0,
        3.0, 4.0
    });

    auto m2 = m1;

    auto m3 = m1 * m2;

    ASSERT_EQ(m3.xx, 7.0);
    ASSERT_EQ(m3.xy, 10.0);
    ASSERT_EQ(m3.yx, 15.0);
    ASSERT_EQ(m3.yy, 22.0);
}

TEST(Linal, MatrixEqual)
{
    gcm::linal::Matrix<2, 3> m1({
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0
    });

    auto m2 = m1;

    auto m3 = m1*2;

    ASSERT_EQ(m1, m2);
    ASSERT_NE(m1, m3);
}

TEST(Linal, MatrixScalarMultiplication)
{
    gcm::linal::Matrix<2, 2, Matrix22Container> m({
        1.0, 2.0,
        3.0, 4.0
    });

    m = 1*m*2;

    ASSERT_EQ(m.xx, 2.0);
    ASSERT_EQ(m.xy, 4.0);
    ASSERT_EQ(m.yx, 6.0);
    ASSERT_EQ(m.yy, 8.0);
}

TEST(Linal, MatrixScalarDivision)
{
    gcm::linal::Matrix<2, 2, Matrix22Container> m({
        2.0, 4.0,
        6.0, 8.0
    });

    m = m/2;

    ASSERT_EQ(m.xx, 1.0);
    ASSERT_EQ(m.xy, 2.0);
    ASSERT_EQ(m.yx, 3.0);
    ASSERT_EQ(m.yy, 4.0);
}

TEST(Linal, MatrixTranspose)
{
    gcm::linal::Matrix<2, 3> m1({
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0
    });

    auto m2 = m1.transpose();

    ASSERT_EQ(m2(0, 0), 1.0);
    ASSERT_EQ(m2(1, 0), 2.0);
    ASSERT_EQ(m2(2, 0), 3.0);
    ASSERT_EQ(m2(0, 1), 4.0);
    ASSERT_EQ(m2(1, 1), 5.0);
    ASSERT_EQ(m2(2, 1), 6.0);

    auto m3 = m2.transpose();

    ASSERT_EQ(m1, m3);
}

TEST(Linal, MatrixTransposeInplace)
{
    gcm::linal::Matrix<2, 2> m1({
        1.0, 2.0,
        3.0, 4.0
    });

    auto m2 = m1;

    m1.transposeInplace();

    ASSERT_EQ(m1(0, 0), 1.0);
    ASSERT_EQ(m1(1, 0), 2.0);
    ASSERT_EQ(m1(0, 1), 3.0);
    ASSERT_EQ(m1(1, 1), 4.0);

    m1.transposeInplace();

    ASSERT_EQ(m1, m2);
}

TEST(Linal, MatrixInvert)
{
    gcm::linal::Matrix<2, 2, Matrix22Container> m1({
        1.0, 2.0,
        1.0, 4.0
    });
    
    gcm::linal::Matrix<2, 2, Matrix22Container> r({
        2.0, -1.0,
        -0.5, 0.5
    });
    
    gcm::linal::Matrix<2, 2> i({
        1.0, 0.0,
        0.0, 1.0
    });


    auto m2 = m1.invert();
    
    ASSERT_EQ(m2 ,r);

    auto m3 = m2.invert();

    ASSERT_EQ(m1, m3);

    ASSERT_EQ(m1*m2, i);
}

TEST(Linal, MatrixInvertInplace)
{
    gcm::linal::Matrix<2, 2, Matrix22Container> m1({
        1.0, 2.0,
        1.0, 4.0
    });

    auto m2 = m1;
    
    gcm::linal::Matrix<2, 2, Matrix22Container> r({
        2.0, -1.0,
        -0.5, 0.5
    });
    
    gcm::linal::Matrix<2, 2> i({
        1.0, 0.0,
        0.0, 1.0
    });


    m1.invertInplace();
    
    ASSERT_EQ(m1 ,r);
    
    ASSERT_EQ(m1*m2, i);

    m1.invertInplace();

    ASSERT_EQ(m1, m2);
}

TEST(Linal, Matrix33Construct)
{
    Matrix33 m1({
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0
    });

    ASSERT_EQ(m1.xx, 1.0);
    ASSERT_EQ(m1.xy, 2.0);
    ASSERT_EQ(m1.xz, 3.0);
    ASSERT_EQ(m1.yx, 4.0);
    ASSERT_EQ(m1.yy, 5.0);
    ASSERT_EQ(m1.yz, 6.0);
    ASSERT_EQ(m1.zx, 7.0);
    ASSERT_EQ(m1.zy, 8.0);
    ASSERT_EQ(m1.zz, 9.0);
}

TEST(Linal, Matrix33Assign)
{
    Matrix33 m1({
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0
    });

    auto m2 = m1;

    ASSERT_EQ(m2.xx, 1.0);
    ASSERT_EQ(m2.xy, 2.0);
    ASSERT_EQ(m2.xz, 3.0);
    ASSERT_EQ(m2.yx, 4.0);
    ASSERT_EQ(m2.yy, 5.0);
    ASSERT_EQ(m2.yz, 6.0);
    ASSERT_EQ(m2.zx, 7.0);
    ASSERT_EQ(m2.zy, 8.0);
    ASSERT_EQ(m2.zz, 9.0);
}

TEST(Linal, Matrix33Determinant)
{
    Matrix33 m({
        1.0, 2.0, 3.0,
        3.0, 2.0, 1.0,
        5.0, 7.0, 11.0
    });

    ASSERT_EQ(determinant(m), -8.0);
}

TEST(Linal, Vector3Construct)
{
    Vector3 m1({
        1.0, 2.0, 3.0,
    });

    ASSERT_EQ(m1.x, 1.0);
    ASSERT_EQ(m1.y, 2.0);
    ASSERT_EQ(m1.z, 3.0);
}

TEST(Linal, Vector3Assign)
{
    Vector3 m1({
        1.0, 2.0, 3.0,
    });

    auto m2 = m1;

    ASSERT_EQ(m2.x, 1.0);
    ASSERT_EQ(m2.y, 2.0);
    ASSERT_EQ(m2.z, 3.0);
}

TEST(Linal, VectorLength)
{
    Vector3 v({
        0.0, 3.0, 4.0,
    });

    ASSERT_EQ(length(v), 5.0);
}

TEST(Linal, VectorDotProduct)
{
    Vector3 v1({
        1.0, 2.0, 3.0
    });

    Vector3 v2({
        -2.0, 4.0, -2.0
    });

    ASSERT_EQ(dotProduct(v1, v2), 0.0);
}

TEST(Linal, VectorCrossProduct)
{
    Vector3 v1({
        1.0, 0.0, 0.0
    });

    Vector3 v2({
        0.0, 1.0, 0.0
    });
    
    Vector3 v3({
        0.0, 0.0, 1.0
    });
    
    Vector3 z({
        0.0, 0.0, 0.0
    });

    ASSERT_EQ(crossProduct(v1, v2), v3);
    ASSERT_EQ(crossProduct(v2, v3), v1);
    ASSERT_EQ(crossProduct(v3, v1), v2);
    
    ASSERT_EQ(crossProduct(v2, v1), -v3);
    ASSERT_EQ(crossProduct(v3, v2), -v1);
    ASSERT_EQ(crossProduct(v1, v3), -v2);
    
    ASSERT_EQ(crossProduct(v1, v1), z);
    ASSERT_EQ(crossProduct(v2, v2), z);
    ASSERT_EQ(crossProduct(v3, v3), z);
}

TEST(Linal, VectorNormalize)
{
    Vector3 v1({
        2.0, 0.0, 0.0
    });

    auto v = normalize(v1);

    ASSERT_NEAR(v.x, 1.0, 1e-5);
    ASSERT_EQ(v.y, 0.0);
    ASSERT_EQ(v.z, 0.0);

    Vector3 v2({
        1.0, 2.0, 3.0        
    });

    v = normalize(v2);

    ASSERT_NEAR(length(v), 1.0, 1e-5);

#if CONFIG_ENABLE_ASSERTIONS
    Vector3 v3({
        0.0, 0.0, 0.0
    });

    ASSERT_THROW(
        normalize(v3),
        Exception
    );
#endif
}

TEST(Linal, RotationMatrix)
{
    Vector3 x({
        1.0, 0.0, 0.0
    });
    
    Vector3 y({
        0.0, 1.0, 0.0
    });
    
    Vector3 z({
        0.0, 0.0, 1.0
    });

    ASSERT_EQ(getXRotationMatrix(M_PI/2)*y, -z);
    ASSERT_EQ(getXRotationMatrix(M_PI/2)*z, y);
    ASSERT_EQ(getXRotationMatrix(-M_PI/2)*y, z);
    ASSERT_EQ(getXRotationMatrix(-M_PI/2)*z, -y);
    
    ASSERT_EQ(getYRotationMatrix(M_PI/2)*x, z);
    ASSERT_EQ(getYRotationMatrix(M_PI/2)*z, -x);
    ASSERT_EQ(getYRotationMatrix(-M_PI/2)*z, x);
    ASSERT_EQ(getYRotationMatrix(-M_PI/2)*x, -z);
    
    ASSERT_EQ(getZRotationMatrix(M_PI/2)*y, x);
    ASSERT_EQ(getZRotationMatrix(M_PI/2)*x, -y);
    ASSERT_EQ(getZRotationMatrix(-M_PI/2)*y, -x);
    ASSERT_EQ(getZRotationMatrix(-M_PI/2)*x, y);
}

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
    Matrix<5, 5> A;
    A(0, 0) = 12; A(1, 1) = 56.333; A(2, 2) = 1; A(3, 3) = 0; A(4, 4) = -34.0022;
    ASSERT_NEAR(A.trace(), 35.3308, EQUALITY_TOLERANCE);
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
    srand((unsigned int)time(0));

    for (int l = 0; l < LINAL_TEST_NUMBER_ITERATIONS; l++) {

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

TEST(Linal, VectorInt)
{
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

TEST(Linal, plainDivision)
{
    VectorInt<3> i = {1, 3, 4};
    Vector<3> r = {0.1, 0.3, 0.7};

    Vector<3> q = plainDivision(r, i);
    ASSERT_NEAR(0.1, q(0), EQUALITY_TOLERANCE);
    ASSERT_NEAR(0.1, q(1), EQUALITY_TOLERANCE);
    ASSERT_NEAR(0.7 / 4, q(2), EQUALITY_TOLERANCE);

    Vector<1> a = {1};
    Vector<1> b = {0};
    ASSERT_NEAR(  std::numeric_limits<real>::max(), plainDivision(  a, b)(0), EQUALITY_TOLERANCE);
    ASSERT_NEAR(- std::numeric_limits<real>::max(), plainDivision(- a, b)(0), EQUALITY_TOLERANCE);
    ASSERT_ANY_THROW(plainDivision(b, b));
}