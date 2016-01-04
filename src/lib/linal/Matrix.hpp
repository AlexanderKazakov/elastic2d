#ifndef LIBGCM_LINAL_MATRIX_HPP
#define LIBGCM_LINAL_MATRIX_HPP

#include <initializer_list>
#include <iostream>
#include <cmath>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


#include "lib/util/Types.hpp"
#include "lib/util/Assertion.hpp"


namespace gcm
{
    namespace linal
    {
        /**
         * Default implementation for matrix values container.
         * We differ SIZE from TM * TN, because SIZE can be bigger in some cases.
         *
         * @tparam TM First matrix dimension (number of strings).
         * @tparam TN Second matrix dimension (number of columns).
         */
        template<int TM, int TN>
        class DefaultMatrixContainer {
        public:
            static const int SIZE = TM * TN; /// size of storage in units of gcm::real
            real values[SIZE];
        };

        struct Matrix22Container {
            static const int SIZE = 2 * 2; /// size of storage in units of gcm::real
            union {
                real values[SIZE];
                struct {
                    real xx;
                    real xy;
                    real yx;
                    real yy;
                };
            };
        };

        template<int TM, typename Container=DefaultMatrixContainer<TM, 1>> class Vector;

        template<class TContainer>
        void clear(TContainer& container) {
            memset(container.values, 0, TContainer::SIZE * sizeof(real));
        };


        /**
         * Generic matrix class.
         *
         * @tparam TM First matrix dimension (number of strings).
         * @tparam TN Second matrix dimension (number of columns).
         * @tparam Container Container class to hold values.
         */
        template<int TM, int TN, typename Container=DefaultMatrixContainer<TM, TN>>
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
            static const int M = TM; /// number of strings
            static const int N = TN; /// number of columns
            /** Default constructor. */
            Matrix();

            /**
             * Copy constructor
             *
             * @param m Matrix to construct from.
             */
            Matrix(const Matrix<TM, TN, Container>& m);

            /**
             * Assignment operator.
             *
             * @param m Matrix to assign values from.
             *
             * @return Reference to modified matrix instance.
             */
            Matrix<TM, TN, Container>& operator=(const Matrix<TM, TN, Container>& m);

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
            Matrix<TN, TM, Container> transpose() const;

            /** Transposes square matrix (modifies matrix contents). */
            void transposeInplace();

            /**
             * Inverses matrix. Note this method returns Matrix<TN, TM> (not Matrix<TM, TN>) to make compilation fail
             * if matrix is not square.
             *
             * @return Inversed matrix.
             */
            Matrix<TN, TM, Container> invert() const;

            /** Inverses matrix modifying its contents. */
            void invertInplace();

            /** @param list of diagonal values. */
            void createDiagonal(const std::initializer_list <real> &list);

            /** Fill in the i-th column. */
            void setColumn(const int i, const Vector<TM> &column);

            /** @return i-th column. */
            Vector<TM> getColumn(const int i) const;

            /** @return values from diagonal multiplied by c in vector */
            Vector<TN> getDiagonalMultipliedBy(const real &c) const;

            /** @return diagonal of this multiplied by B in vector */
            Vector<TM> diagonalMultiply(const Matrix<TN, TM> &B) const;

            /** @return trace of the matrix */
            real trace() const;
        };

        template<int TM, int TN, typename Container>
        Matrix<TM, TN, Container>::Matrix()
        {
            static_assert(sizeof(this->values) >= sizeof(real)*TM*TN, "Container must have enough memory to store values");
        }

        template<int TM, int TN, typename Container>
        Matrix<TM, TN, Container>::Matrix(std::initializer_list<real> values): Matrix()
        {
            int i = 0;
            for (auto value: values)
                this->values[i++] = value;
        }
        
        template<int TM, int TN, typename Container>
        Matrix<TM, TN, Container>::Matrix(const Matrix<TM, TN, Container>& m)
        {
            (*this) = m;
        }

        template<int TM, int TN, typename Container>
        Matrix<TM, TN, Container>& Matrix<TM, TN, Container>::operator=(const Matrix<TM, TN, Container>& m)
        {
            for (int i = 0; i < TM; i++)
                for (int j = 0; j < TN; j++)
                    (*this)(i, j) = m(i, j);

            return *this;
        }

        template<int TM, int TN, typename Container>
        inline
        int Matrix<TM, TN, Container>::getIndex(int i, int j) const
        {
            assert_lt(i, TM);
            assert_lt(j, TN);

            return i*TN+j;
        }

        template<int TM, int TN, typename Container>
        inline
        real Matrix<TM, TN, Container>::operator()(const int i, const int j) const
        {
            return this->values[getIndex(i, j)];
        }

        template<int TM, int TN, typename Container>
        inline
        real& Matrix<TM, TN, Container>::operator()(const int i, const int j)
        {
            return this->values[getIndex(i, j)];
        }

        template<int TM, int TN, typename Container>
        Matrix<TN, TM, Container> Matrix<TM, TN, Container>::transpose() const
        {
            Matrix<TN, TM, Container> result;

            for (int i = 0; i < TM; i++)
                for (int j = 0; j < TN; j++)
                    result(j, i) = (*this)(i, j);

            return result;
        }
        
        template<int TM, int TN, typename Container>
        void Matrix<TM, TN, Container>::transposeInplace()
        {
            (*this) = transpose();
        }
        
        template<int TM, int TN, typename Container>
        Matrix<TN, TM, Container> Matrix<TM, TN, Container>::invert() const
        {
            Matrix<TN, TM, Container> result;

            gsl_matrix* Z1 = gsl_matrix_alloc(TM, TM);
            gsl_matrix* Z = gsl_matrix_alloc(TM, TM);
            gsl_permutation* perm = gsl_permutation_alloc(TM);
            int k;

            for (int i = 0; i < TM; i++)
                for (int j = 0; j < TM; j++)
                    gsl_matrix_set(Z1, i, j, (*this)(i, j));

            if (gsl_linalg_LU_decomp(Z1, perm, &k))
                THROW_INVALID_ARG("gsl_linalg_LU_decomp failed");
            
            if (gsl_linalg_LU_invert(Z1, perm, Z))
                THROW_INVALID_ARG("gsl_linalg_LU_invert failed");
            
            for (int i = 0; i < TM; i++)
                for (int j = 0; j < TM; j++)
                    result(i, j) = gsl_matrix_get(Z, i, j);

            gsl_permutation_free(perm);
            gsl_matrix_free(Z);
            gsl_matrix_free(Z1);

            return result;
        }
        
        template<int TM, int TN, typename Container>
        void Matrix<TM, TN, Container>::invertInplace()
        {
            (*this) = invert();
        }

        template<int TM, int TN, typename Container>
        void Matrix<TM, TN, Container>::createDiagonal(const std::initializer_list<real>& list) {
            clear(*this);
            int i = 0;
            for (auto& value : list) {
                (*this)(i, i++) = value;
            }
        }

        template<int TM, int TN, typename Container>
        void Matrix<TM, TN, Container>::setColumn(const int i, const Vector<TM>& column) {
            for (int j = 0; j < TM; j++) {
                (*this)(j, i) = column(j);
            }
        }

        template<int TM, int TN, typename Container>
        Vector<TM> Matrix<TM, TN, Container>::getColumn(const int i) const {
            Vector<TM> ans;
            for (int j = 0; j < TM; j++) {
                ans(j) = (*this)(j, i);
            }
            return ans;
        }

        template<int TM, int TN, typename Container>
        Vector<TN> Matrix<TM, TN, Container>::getDiagonalMultipliedBy(const real &c) const {
            Vector<TN> ans;
            for (int i = 0; i < TN; i++) {
                ans(i) = (*this)(i, i) * c;
            }
            return ans;
        }

        template<int TM, int TN, typename Container>
        Vector<TM> Matrix<TM, TN, Container>::diagonalMultiply(const Matrix<TN, TM> &B) const {
            Vector<TM> ans;
            for (int i = 0; i < TM; i++) {
                ans(i) = 0;
                for (int j = 0; j < TM; j++) {
                    ans(i) += (*this)(i, j) * B(j, i);
                }
            }
            return ans;
        }

        /*template<int TM, int TN, typename Container>
        Vector<TM> Matrix<TM, TN, Container>::operator*(const Vector<TN> &b) const {
            Vector<TN> c; // c = this * b
            for (int i = 0; i < TN; i++) {
                c(i) = 0.0;
                for (int j = 0; j < TM; j++) {
                    c(i) += (*this)(i, j) * b(j);
                }
            }
            return c;
        }*/

        template<int TM, int TN, typename Container>
        real Matrix<TM, TN, Container>::trace() const {
            real trace = 0.0;
            for (int i = 0; i < TN; i++) {
                trace += (*this)(i, i);
            }
            return trace;
        }

        /**
         * Computes negative of matrix B
         *
         * @tparam TM First matrix dimension
         * @tparam TN Second matrix dimension
         * @tparam Container Matrix container type
         * @param m Matrix to compute negative for.
         *
         * @return Negative of matrix
         */
        template<int TM, int TN, typename Container>
        Matrix<TM, TN, Container> operator-(const Matrix<TM, TN, Container>& m)
        {
            Matrix<TM, TN, Container> result;

            for (int i = 0; i < TM; i++)
                for (int j = 0; j < TN; j++)
                    result(i, j) = -m(i, j);

            return result;
        }

        /**
         * Computes summ of two matrixes. Generic implementation.
         *
         * @tparam TM First matrix dimension.
         * @tparam TN Second matrix dimension.
         * @tparam Container1 Container type for first summ item.
         * @tparam Container2 Container type for second  summ item.
         * @tparam Container3 Container type for result.
         * @param m1 First summ item.
         * @param m2 Second summ item.
         *
         * @return Matrix summ (m1+m2).
         */
        template<int TM, int TN, typename Container1, typename Container2, typename Container3>
        Matrix<TM, TN, Container3> operator+(const Matrix<TM, TN, Container1>& m1, const Matrix<TM, TN, Container2>& m2)
        {
            Matrix<TM, TN, Container3> result;

            for (int i = 0; i < TM; i++)
                for (int j = 0; j < TN; j++)
                    result(i, j) = m1(i, j) + m2(i, j);

            return result;
        }
        
        /**
         * Computes summ of two matrixes. Most useful specialization, made in assumption that result matrix should
         * use default matrix container.
         *
         * @tparam TM First matrix dimension.
         * @tparam TN Second matrix dimension.
         * @tparam Container1 Container type for first summ item.
         * @tparam Container2 Container type for second  summ item.
         * @param m1 First summ item.
         * @param m2 Second summ item.
         *
         * @return Matrix summ (m1+m2).
         */
        template<int TM, int TN, typename Container1, typename Container2>
        Matrix<TM, TN, DefaultMatrixContainer<TM, TN>> operator+(const Matrix<TM, TN, Container1>& m1, const Matrix<TM, TN, Container2>& m2)
        {
            return operator+<TM, TN, Container1, Container2, DefaultMatrixContainer<TM, TN>>(m1, m2);
        }

        /**
         * Computes summ of two matrixes. Specialized implementation, made in assumption that all matrixes
         * use default matrix container.
         *
         * @tparam TM First matrix dimension.
         * @tparam TN Second matrix dimension.
         * @tparam Container Container type for all matrixes.
         * @param m1 First summ item.
         * @param m2 Second summ item.
         *
         * @return Matrix summ (m1+m2).
         */
        template<int TM, int TN, typename Container>
        Matrix<TM, TN, Container> operator+(const Matrix<TM, TN, Container>& m1, const Matrix<TM, TN, Container>& m2)
        {
            return operator+<TM, TN, Container, Container, Container>(m1, m2);
        }
        
        /**
         * Computes difference of two matrixes. Generic implementation.
         *
         * @tparam TM First matrix dimension.
         * @tparam TN Second matrix dimension.
         * @tparam Container1 Container type for first difference item.
         * @tparam Container2 Container type for second  difference item.
         * @tparam Container3 Container type for result.
         * @param m1 First difference item.
         * @param m2 Second difference item.
         *
         * @return Matrix difference (m1-m2).
         */
        template<int TM, int TN, typename Container1, typename Container2, typename Container3>
        Matrix<TM, TN, Container3> operator-(const Matrix<TM, TN, Container1>& m1, const Matrix<TM, TN, Container2>& m2)
        {
            Matrix<TM, TN, Container3> result;

            for (int i = 0; i < TM; i++)
                for (int j = 0; j < TN; j++)
                    result(i, j) = m1(i, j) - m2(i, j);

            return result;
        }
        
        /**
         * Computes difference of two matrixes. Most useful specialization, made in assumption that result matrix should
         * use default matrix container.
         *
         * @tparam TM First matrix dimension.
         * @tparam TN Second matrix dimension.
         * @tparam Container1 Container type for first difference item.
         * @tparam Container2 Container type for second  difference item.
         * @param m1 First difference item.
         * @param m2 Second difference item.
         *
         * @return Matrix difference (m1-m2).
         */
        template<int TM, int TN, typename Container1, typename Container2>
        Matrix<TM, TN, DefaultMatrixContainer<TM, TN>> operator-(const Matrix<TM, TN, Container1>& m1, const Matrix<TM, TN, Container2>& m2)
        {
            return operator-<TM, TN, Container1, Container2, DefaultMatrixContainer<TM, TN>>(m1, m2);
        }

        /**
         * Computes difference of two matrixes. Specialized implementation, made in assumption that all matrixes
         * use default matrix container.
         *
         * @tparam TM First matrix dimension.
         * @tparam TN Second matrix dimension.
         * @tparam Container Container type for all matrixes.
         * @param m1 First difference item.
         * @param m2 Second difference item.
         *
         * @return Matrix difference (m1-m2).
         */
        template<int TM, int TN, typename Container>
        Matrix<TM, TN, Container> operator-(const Matrix<TM, TN, Container>& m1, const Matrix<TM, TN, Container>& m2)
        {
            return operator-<TM, TN, Container, Container, Container>(m1, m2);
        }
        
        /**
         * Computes product of two matrixes. Generic implementation.
         *
         * @tparam TM First dimension of first matrix (TM x TN).
         * @tparam TN First (second) dimension of second (first) matrix (TM x TN or TN x TK respectively).
         * @tparam TK Second dimension of second matrix (TN x TK).
         * @tparam Container1 Container type for first product item.
         * @tparam Container2 Container type for second  product item.
         * @tparam Container3 Container type for result.
         * @param m1 First product item.
         * @param m2 Second product item.
         *
         * @return Matrix product (m1*m2).
         */
        template<int TM, int TN, int TK, typename Container1, typename Container2, typename Container3>
        Matrix<TM, TK, Container3> operator*(const Matrix<TM, TN, Container1>& m1, const Matrix<TN, TK, Container2>& m2)
        {
            Matrix<TM, TK, Container3> result;

            for (int i = 0; i < TM; i++)
                
                for (int j = 0; j < TK; j++)
                {
                    result(i, j) = 0;
                    for (int n = 0; n < TN; n++)
                        result(i, j) += m1(i, n) * m2(n, j);
                }

            return result;
        }
        
        /**
         * Computes product of two matrixes. Most useful specialization, made in assumption that result matrix should
         * use default matrix container.
         *
         * @tparam TM First dimension of first matrix (TM x TN).
         * @tparam TN First (second) dimension of second (first) matrix (TM x TN or TN x TK respectively).
         * @tparam TK Second dimension of second matrix (TN x TK).
         * @tparam Container1 Container type for first product item.
         * @tparam Container2 Container type for second  product item.
         * @param m1 First product item.
         * @param m2 Second product item.
         *
         * @return Matrix product (m1*m2).
         */
        template<int TM, int TN, int TK, typename Container1, typename Container2>
        Matrix<TM, TK, DefaultMatrixContainer<TM, TK>> operator*(const Matrix<TM, TN, Container1>& m1, const Matrix<TN, TK, Container2>& m2)
        {
            return operator*<TM, TN, TK, Container1, Container2, DefaultMatrixContainer<TM, TK>>(m1, m2);
        }
        
        /**
         * Computes product of two matrixes. Most useful specialization, made in assumption that both matrixes are
         * square and have the same container type.
         *
         * @tparam TM Matrix size.
         * @tparam Container Container type for first product item.
         * @param m1 First product item.
         * @param m2 Second product item.
         *
         * @return Matrix product (m1*m2).
         */
        template<int TM, typename Container>
        Matrix<TM, TM, Container> operator*(const Matrix<TM, TM, Container>& m1, const Matrix<TM, TM, Container>& m2)
        {
            return operator*<TM, TM, TM, Container, Container, Container>(m1, m2);
        }

        /**
         * Performs scalar multiplication.
         *
         * @tparam TM First matrix dimesion.
         * @tparam TN Second matrix dimension.
         * @tparam Container Container type of matrix.
         * @param m Matrix to multiply.
         * @param x Scalar to multiply by.
         *
         * @return Result of scalar multiplication.
         */
        template<int TM, int TN, typename Container>
        Matrix<TM, TN, Container> operator*(const Matrix<TM, TN, Container>& m, const real x)
        {
            Matrix<TM, TN, Container> result;

            for (int i = 0; i < TM; i++)
                for (int j = 0; j < TN; j++)
                    result(i, j) = m(i, j) * x;

            return result;
        }
        
        /**
         * Performs scalar multiplication.
         *
         * @tparam TM First matrix dimesion.
         * @tparam TN Second matrix dimension.
         * @tparam Container Container type of matrix.
         * @param x Scalar to multiply by.
         * @param m Matrix to multiply.
         *
         * @return Result of scalar multiplication.
         */
        template<int TM, int TN, typename Container>
        Matrix<TM, TN, Container> operator*(const real x, const Matrix<TM, TN, Container>& m)
        {
            return m * x;
        }
        
        /**
         * Performs scalar division.
         *
         * @tparam TM First matrix dimesion.
         * @tparam TN Second matrix dimension.
         * @tparam Container Container type of matrix.
         * @param m Matrix to divide.
         * @param x Scalar to divide by.
         *
         * @return Result of scalar division.
         */
        template<int TM, int TN, typename Container>
        Matrix<TM, TN, Container> operator/(const Matrix<TM, TN, Container>& m, const real x)
        {
            return m * (1/x);
        }

        /**
         * Add m2 to m1. Generic implementation.
         *
         * @tparam TM First matrix dimension.
         * @tparam TN Second matrix dimension.
         * @tparam Container1 Container type for first summ item.
         * @tparam Container2 Container type for second  summ item.
         * @param m1 First summ item.
         * @param m2 Second summ item.
         */
        template<int TM, int TN, typename Container1, typename Container2>
        void operator+=(Matrix<TM, TN, Container1>& m1, const Matrix<TM, TN, Container2>& m2)
        {
            for (int i = 0; i < TM; i++)
                for (int j = 0; j < TN; j++)
                    m1(i, j) = m1(i, j) + m2(i, j);
        }

        template<int TM, int TN, typename Container1, typename Container2>
        bool operator==(const Matrix<TM, TN, Container1>& m1, const Matrix<TM, TN, Container2>& m2)
        {
            for (int i = 0; i < TM; i++)
                for (int j = 0; j < TN; j++)
                    // FIXME Should this constant be replaced by something context-specific?
                    if (fabs(m1(i, j) - m2(i, j)) > EQUALITY_TOLERANCE)
                        return false;

            return true;
        }
        
        template<int TM, int TN, typename Container1, typename Container2>
        bool operator!=(const Matrix<TM, TN, Container1>& m1, const Matrix<TM, TN, Container2>& m2)
        {
            return !(m1 == m2);
        }
    };
};

namespace std {
    template<int TM, int TN>
    inline std::ostream& operator<<(std::ostream &os, const gcm::linal::Matrix<TM, TN>& matrix) {

        os << "Matrix:\n";
        for (int i = 0; i < TM; i++) {
            for (int j = 0; j < TN; j++) {
                os << matrix(i, j) << "\t";
            }
            os << "\n";
        }

        return os;
    };
};

#endif // LIBGCM_LINAL_MATRIX_HPP
