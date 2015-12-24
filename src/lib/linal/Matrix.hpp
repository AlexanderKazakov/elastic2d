#ifndef LIBGCM_LINAL_MATRIX_HPP
#define LIBGCM_LINAL_MATRIX_HPP

#include "lib/util/Types.hpp"
#include "lib/util/Assertion.hpp"

#include <initializer_list>
#include <cmath>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


namespace gcm
{
    namespace linal
    {
        /**
         * Default implementation for matrix values container.
         *
         * @tparam M First matrix dimension.
         * @tparam N Second matrix dimension.
         */
        template<int M, int N>
        class DefaultMatrixContainer
        {
            public:
                real values[M*N];
        };

        struct Matrix22Container
        {
            union
            {
                real values[4];
                struct
                {
                    real xx;
                    real xy;
                    real yx;
                    real yy;
                };
            };
        };

        /**
         * Generic matrix class.
         *
         * @tparam M First matrix dimension.
         * @tparam N Second matrix dimension.
         * @tparam Container Container class to hold values.
         */
        template<int M, int N, typename Container=DefaultMatrixContainer<M, N>>
        class Matrix: public Container
        {
            protected:
                /**
                 * Returns values array index of matrix component.
                 *
                 * @param i First index.
                 * @param j Second index.
                 *
                 * @return Values array index.
                 */
                int getIndex(const int i, const int j) const;
            public:
                /**
                 * Default constructor.
                 */
                Matrix();

                /**
                 * Copy constructor
                 *
                 * @param m Matrix to construct from.
                 */
                Matrix(const Matrix<M, N, Container>& m);

                /**
                 * Assignment operator.
                 *
                 * @param m Matrix to assign values from.
                 *
                 * @return Reference to modified matrix instance.
                 */
                Matrix<M, N, Container>& operator=(const Matrix<M, N, Container>& m);

                /**
                 * Constructor that initializes matrix with specified values.
                 *
                 * @param values Values to initialize matrix with.
                 */
                Matrix(std::initializer_list<real> values);
                
                /**
                 * Returns matrix component.
                 *
                 * @param i First index.
                 * @param j Second index.
                 *
                 * @return Corresponding matrix component.
                 */
                real operator()(const int i, const int j) const;
                
                /**
                 * Returns reference to matrix component, used to modify matrix content.
                 *
                 * @param i First index.
                 * @param j Second index.
                 *
                 * @return Reference to corresponding matrix component.
                 */
                real& operator()(const int i, const int j);

                /**
                 * Transposes matrix.
                 *
                 * @return Transposed matrix.
                 */
                Matrix<N, M, Container> transpose() const;

                /**
                 * Transposes square matrix (modifies matrix contents).
                 */
                void transposeInplace();

                /**
                 * Inverses matrix. Note this method returns Matrix<N, M> (not Matrix<M, N>) to make compilation fail
                 * if matrix is not square.
                 *
                 * @return Inversed matrix.
                 */
                Matrix<N, M, Container> invert() const;

                /**
                 * Inverses matrix modifying its contents.
                 */
                void invertInplace();
        };

        template<int M, int N, typename Container>
        Matrix<M, N, Container>::Matrix()
        {
            static_assert(sizeof(this->values) >= sizeof(real)*M*N, "Container must have enough memory to store values");
        }

        template<int M, int N, typename Container>
        Matrix<M, N, Container>::Matrix(std::initializer_list<real> values): Matrix()
        {
            int i = 0;
            for (auto value: values)
                this->values[i++] = value;
        }
        
        template<int M, int N, typename Container>
        Matrix<M, N, Container>::Matrix(const Matrix<M, N, Container>& m)
        {
            (*this) = m;
        }

        template<int M, int N, typename Container>
        Matrix<M, N, Container>& Matrix<M, N, Container>::operator=(const Matrix<M, N, Container>& m)
        {
            for (int i = 0; i < M; i++)
                for (int j = 0; j < N; j++)
                    (*this)(i, j) = m(i, j);

            return *this;
        }

        template<int M, int N, typename Container>
        inline
        int Matrix<M, N, Container>::getIndex(int i, int j) const
        {
            assert_lt(i, M);
            assert_lt(j, N);

            return i*N+j;
        }

        template<int M, int N, typename Container>
        inline
        gcm::real Matrix<M, N, Container>::operator()(const int i, const int j) const
        {
            return this->values[getIndex(i, j)];
        }

        template<int M, int N, typename Container>
        inline
        gcm::real& Matrix<M, N, Container>::operator()(const int i, const int j)
        {
            return this->values[getIndex(i, j)];
        }

        template<int M, int N, typename Container>
        Matrix<N, M, Container> Matrix<M, N, Container>::transpose() const
        {
            Matrix<N, M, Container> result;

            for (int i = 0; i < M; i++)
                for (int j = 0; j < N; j++)
                    result(j, i) = (*this)(i, j);

            return result;
        }
        
        template<int M, int N, typename Container>
        void Matrix<M, N, Container>::transposeInplace()
        {
            (*this) = transpose();
        }
        
        template<int M, int N, typename Container>
        Matrix<N, M, Container> Matrix<M, N, Container>::invert() const
        {
            Matrix<N, M, Container> result;

            gsl_matrix* Z1 = gsl_matrix_alloc(M, M);
            gsl_matrix* Z = gsl_matrix_alloc(M, M);
            gsl_permutation* perm = gsl_permutation_alloc(M);
            int k;

            for (int i = 0; i < M; i++)
                for (int j = 0; j < M; j++)
                    gsl_matrix_set(Z1, i, j, (*this)(i, j));

            if (gsl_linalg_LU_decomp(Z1, perm, &k))
                THROW_INVALID_ARG("gsl_linalg_LU_decomp failed");
            
            if (gsl_linalg_LU_invert(Z1, perm, Z))
                THROW_INVALID_ARG("gsl_linalg_LU_invert failed");
            
            for (int i = 0; i < M; i++)
                for (int j = 0; j < M; j++)
                    result(i, j) = gsl_matrix_get(Z, i, j);

            gsl_permutation_free(perm);
            gsl_matrix_free(Z);
            gsl_matrix_free(Z1);

            return result;
        }
        
        template<int M, int N, typename Container>
        void Matrix<M, N, Container>::invertInplace()
        {
            (*this) = invert();
        }
        
        /**
         * Computes negative of matrix B
         *
         * @tparam M First matrix dimension
         * @tparam N Second matrix dimension
         * @tparam Container Matrix container type
         * @param m Matrix to compute negative for.
         *
         * @return Negative of matrix
         */
        template<int M, int N, typename Container>
        Matrix<M, N, Container> operator-(const Matrix<M, N, Container>& m)
        {
            Matrix<M, N, Container> result;

            for (int i = 0; i < M; i++)
                for (int j = 0; j < N; j++)
                    result(i, j) = -m(i, j);

            return result;
        }

        /**
         * Computes summ of two matrixes. Generic implementation.
         *
         * @tparam M First matrix dimension.
         * @tparam N Second matrix dimension.
         * @tparam Container1 Container type for first summ item.
         * @tparam Container2 Container type for second  summ item.
         * @tparam Container3 Container type for result.
         * @param m1 First summ item.
         * @param m2 Second summ item.
         *
         * @return Matrix summ (m1+m2).
         */
        template<int M, int N, typename Container1, typename Container2, typename Container3>
        Matrix<M, N, Container3> operator+(const Matrix<M, N, Container1>& m1, const Matrix<M, N, Container2>& m2)
        {
            Matrix<M, N, Container3> result;

            for (int i = 0; i < M; i++)
                for (int j = 0; j < N; j++)
                    result(i, j) = m1(i, j) + m2(i, j);

            return result;
        }
        
        /**
         * Computes summ of two matrixes. Most useful specialization, made in assumption that result matrix should
         * use default matrix container.
         *
         * @tparam M First matrix dimension.
         * @tparam N Second matrix dimension.
         * @tparam Container1 Container type for first summ item.
         * @tparam Container2 Container type for second  summ item.
         * @param m1 First summ item.
         * @param m2 Second summ item.
         *
         * @return Matrix summ (m1+m2).
         */
        template<int M, int N, typename Container1, typename Container2>
        Matrix<M, N, DefaultMatrixContainer<M, N>> operator+(const Matrix<M, N, Container1>& m1, const Matrix<M, N, Container2>& m2)
        {
            return operator+<M, N, Container1, Container2, DefaultMatrixContainer<M, N>>(m1, m2);
        }

        /**
         * Computes summ of two matrixes. Specialized implementation, made in assumption that all matrixes
         * use default matrix container.
         *
         * @tparam M First matrix dimension.
         * @tparam N Second matrix dimension.
         * @tparam Container Container type for all matrixes.
         * @param m1 First summ item.
         * @param m2 Second summ item.
         *
         * @return Matrix summ (m1+m2).
         */
        template<int M, int N, typename Container>
        Matrix<M, N, Container> operator+(const Matrix<M, N, Container>& m1, const Matrix<M, N, Container>& m2)
        {
            return operator+<M, N, Container, Container, Container>(m1, m2);
        }
        
        /**
         * Computes difference of two matrixes. Generic implementation.
         *
         * @tparam M First matrix dimension.
         * @tparam N Second matrix dimension.
         * @tparam Container1 Container type for first difference item.
         * @tparam Container2 Container type for second  difference item.
         * @tparam Container3 Container type for result.
         * @param m1 First difference item.
         * @param m2 Second difference item.
         *
         * @return Matrix difference (m1-m2).
         */
        template<int M, int N, typename Container1, typename Container2, typename Container3>
        Matrix<M, N, Container3> operator-(const Matrix<M, N, Container1>& m1, const Matrix<M, N, Container2>& m2)
        {
            Matrix<M, N, Container3> result;

            for (int i = 0; i < M; i++)
                for (int j = 0; j < N; j++)
                    result(i, j) = m1(i, j) - m2(i, j);

            return result;
        }
        
        /**
         * Computes difference of two matrixes. Most useful specialization, made in assumption that result matrix should
         * use default matrix container.
         *
         * @tparam M First matrix dimension.
         * @tparam N Second matrix dimension.
         * @tparam Container1 Container type for first difference item.
         * @tparam Container2 Container type for second  difference item.
         * @param m1 First difference item.
         * @param m2 Second difference item.
         *
         * @return Matrix difference (m1-m2).
         */
        template<int M, int N, typename Container1, typename Container2>
        Matrix<M, N, DefaultMatrixContainer<M, N>> operator-(const Matrix<M, N, Container1>& m1, const Matrix<M, N, Container2>& m2)
        {
            return operator-<M, N, Container1, Container2, DefaultMatrixContainer<M, N>>(m1, m2);
        }

        /**
         * Computes difference of two matrixes. Specialized implementation, made in assumption that all matrixes
         * use default matrix container.
         *
         * @tparam M First matrix dimension.
         * @tparam N Second matrix dimension.
         * @tparam Container Container type for all matrixes.
         * @param m1 First difference item.
         * @param m2 Second difference item.
         *
         * @return Matrix difference (m1-m2).
         */
        template<int M, int N, typename Container>
        Matrix<M, N, Container> operator-(const Matrix<M, N, Container>& m1, const Matrix<M, N, Container>& m2)
        {
            return operator-<M, N, Container, Container, Container>(m1, m2);
        }
        
        /**
         * Computes product of two matrixes. Generic implementation.
         *
         * @tparam M First dimension of first matrix (M x N).
         * @tparam N First (second) dimension of second (first) matrix (M x N or N x K respectively).
         * @tparam K Second dimension of second matrix (N x K).
         * @tparam Container1 Container type for first product item.
         * @tparam Container2 Container type for second  product item.
         * @tparam Container3 Container type for result.
         * @param m1 First product item.
         * @param m2 Second product item.
         *
         * @return Matrix product (m1*m2).
         */
        template<int M, int N, int K, typename Container1, typename Container2, typename Container3>
        Matrix<M, K, Container3> operator*(const Matrix<M, N, Container1>& m1, const Matrix<N, K, Container2>& m2)
        {
            Matrix<M, K, Container3> result;

            for (int i = 0; i < M; i++)
                
                for (int j = 0; j < K; j++)
                {
                    result(i, j) = 0;
                    for (int n = 0; n < N; n++)
                        result(i, j) += m1(i, n) * m2(n, j);
                }

            return result;
        }
        
        /**
         * Computes product of two matrixes. Most useful specialization, made in assumption that result matrix should
         * use default matrix container.
         *
         * @tparam M First dimension of first matrix (M x N).
         * @tparam N First (second) dimension of second (first) matrix (M x N or N x K respectively).
         * @tparam K Second dimension of second matrix (N x K).
         * @tparam Container1 Container type for first product item.
         * @tparam Container2 Container type for second  product item.
         * @param m1 First product item.
         * @param m2 Second product item.
         *
         * @return Matrix product (m1*m2).
         */
        template<int M, int N, int K, typename Container1, typename Container2>
        Matrix<M, K, DefaultMatrixContainer<M, K>> operator*(const Matrix<M, N, Container1>& m1, const Matrix<N, K, Container2>& m2)
        {
            return operator*<M, N, K, Container1, Container2, DefaultMatrixContainer<M, K>>(m1, m2);
        }
        
        /**
         * Computes product of two matrixes. Most useful specialization, made in assumption that both matrixes are
         * square and have the same container type.
         *
         * @tparam M Matrix size.
         * @tparam Container Container type for first product item.
         * @param m1 First product item.
         * @param m2 Second product item.
         *
         * @return Matrix product (m1*m2).
         */
        template<int M, typename Container>
        Matrix<M, M, Container> operator*(const Matrix<M, M, Container>& m1, const Matrix<M, M, Container>& m2)
        {
            return operator*<M, M, M, Container, Container, Container>(m1, m2);
        }

        /**
         * Performs scalar multiplication.
         *
         * @tparam M First matrix dimesion.
         * @tparam N Second matrix dimension.
         * @tparam Container Container type of matrix.
         * @param m Matrix to multiply.
         * @param x Scalar to multiply by.
         *
         * @return Result of scalar multiplication.
         */
        template<int M, int N, typename Container>
        Matrix<M, N, Container> operator*(const Matrix<M, N, Container>& m, const real x)
        {
            Matrix<M, N, Container> result;

            for (int i = 0; i < M; i++)
                for (int j = 0; j < N; j++)
                    result(i, j) = m(i, j) * x;

            return result;
        }
        
        /**
         * Performs scalar multiplication.
         *
         * @tparam M First matrix dimesion.
         * @tparam N Second matrix dimension.
         * @tparam Container Container type of matrix.
         * @param x Scalar to multiply by.
         * @param m Matrix to multiply.
         *
         * @return Result of scalar multiplication.
         */
        template<int M, int N, typename Container>
        Matrix<M, N, Container> operator*(const real x, const Matrix<M, N, Container>& m)
        {
            return m * x;
        }
        
        /**
         * Performs scalar division.
         *
         * @tparam M First matrix dimesion.
         * @tparam N Second matrix dimension.
         * @tparam Container Container type of matrix.
         * @param m Matrix to divide.
         * @param x Scalar to divide by.
         *
         * @return Result of scalar division.
         */
        template<int M, int N, typename Container>
        Matrix<M, N, Container> operator/(const Matrix<M, N, Container>& m, const real x)
        {
            return m * (1/x);
        }

        template<int M, int N, typename Container1, typename Container2>
        bool operator==(const Matrix<M, N, Container1>& m1, const Matrix<M, N, Container2>& m2)
        {
            for (int i = 0; i < M; i++)
                for (int j = 0; j < N; j++)
                    // FIXME Should this constant be replaced by something context-specific?
                    if (fabs(m1(i, j) - m2(i, j)) > EQUALITY_TOLERANCE)
                        return false;

            return true;
        }
        
        template<int M, int N, typename Container1, typename Container2>
        bool operator!=(const Matrix<M, N, Container1>& m1, const Matrix<M, N, Container2>& m2)
        {
            return !(m1 == m2);
        }
    };
};

#endif // LIBGCM_LINAL_MATRIX_HPP
