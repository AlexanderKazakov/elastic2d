#ifndef LIBGCM_LINAL_MATRIX_HPP
#define LIBGCM_LINAL_MATRIX_HPP

#include <initializer_list>
#include <iostream>
#include <string.h>
#include <cmath>

#include <lib/util/Types.hpp>
#include <lib/util/Assertion.hpp>
#include <lib/util/GslUtils.hpp>


namespace gcm {
namespace linal {
/**
 * Default implementation for matrix values container.
 * @tparam TM First matrix dimension (number of strings).
 * @tparam TN Second matrix dimension (number of columns).
 */
template<int TM, int TN>
class DefaultMatrixContainer {
public:
	static const int SIZE = TM * TN; ///< size of storage in units of gcm::real
	real values[SIZE];
};


/**
 * Generic matrix class.
 *
 * @tparam TM First matrix dimension (number of strings).
 * @tparam TN Second matrix dimension (number of columns).
 * @tparam Container Container class to hold values.
 */
template<int TM, int TN, typename Container = DefaultMatrixContainer<TM, TN> >
class Matrix : public Container {
public:
	typedef Container ContainerType;
	static const int M = TM; ///< number of strings
	static const int N = TN; ///< number of columns

	static Matrix zeros() {
		Matrix m;
		return clear(m);
	}

	Matrix() { }
	Matrix& operator=(const Matrix& m2) = default;

	/** @param values Values to initialize matrix with, string by string */
	Matrix(std::initializer_list<real> list) {
		this->initialize(list);
	}

	/**
	 * Copy constructor
	 * @param m Matrix to construct from
	 */
	template<typename Container2>
	Matrix(const Matrix<TM, TN, Container2>& m2) {
		static_assert(this->SIZE == TM * TN,
		              "Container must have enough memory to store values");
		(*this) = m2;
	}

	/**
	 * Assignment operator from matrix of equal size and any container
	 * @return reference to modified matrix instance
	 */
	template<typename Container2>
	Matrix& operator=(const Matrix<TM, TN, Container2>& m2) {
		// TODO - rvalues?
		static_assert(this->SIZE == m2.SIZE, "Containers must have equal size");
		memcpy(this->values, m2.values, sizeof(this->values));
		return *this;
	}

	/** @param values Values to initialize matrix with, string by string */
	void initialize(std::initializer_list<real> list);

	/** Read-only access to matrix component */
	real operator()(const int i, const int j) const {
		return this->values[getIndex(i, j)];
	}

	/** Read/write access to matrix component */
	real& operator()(const int i, const int j) {
		return this->values[getIndex(i, j)];
	}

	/** Read-only access to vector component */
	real operator()(const int i) const {
		return this->values[i];
	}

	/** Read/write access to vector component */
	real& operator()(const int i) {
		return this->values[i];
	}

	/** @return transposed matrix */
	Matrix<TN, TM, Container> transpose() const;

	/** Transposes square matrix (modifying matrix contents) */
	void transposeInplace();

	/**
	 * Inverses matrix. Note this method returns Matrix<TN, TM> (not Matrix<TM, TN>) to make
	 *compilation fail
	 * if matrix is not square.
	 * @return inverted matrix.
	 */
	Matrix<TN, TM, Container> invert() const;

	/** Inverses matrix modifying its contents */
	void invertInplace();

	/** @return i-th column. */
	template<typename Container2 = DefaultMatrixContainer<TM, 1> >
	Matrix<TM, 1, Container2> getColumn(const int i) const {
		Matrix<TM, 1, Container2> ans;
		for (int j = 0; j < TM; j++) {
			ans(j) = (*this)(j, i);
		}
		return ans;
	}

	/** set i-th column */
	template<typename Container2>
	void setColumn(const int i, const Matrix<TM, 1, Container2>& column) {
		for (int j = 0; j < TM; j++) {
			(*this)(j, i) = column(j);
		}
	}
	
	/** @return i-th string. */
	template<typename Container2 = DefaultMatrixContainer<TN, 1> >
	Matrix<TN, 1, Container2> getString(const int i) const {
		Matrix<TN, 1, Container2> ans;
		for (int j = 0; j < TN; j++) {
			ans(j) = (*this)(i, j);
		}
		return ans;
	}
	
	/** set i-th string */
	template<typename Container2>
	void setString(const int i, const Matrix<TN, 1, Container2>& column) {
		for (int j = 0; j < TN; j++) {
			(*this)(i, j) = column(j);
		}
	}

	/** @return in vector diagonal of this matrix multiplied by matrix B */
	Matrix<TM, 1,
	       DefaultMatrixContainer<TM, 1> > diagonalMultiply(const Matrix<TN, TM>& B) const;

	/** @return trace of the matrix */
	real trace() const;

protected:
	/** @return values array index of matrix component */
	int getIndex(const int i, const int j) const {
		return i * TN + j;
	}

};


template<int TM, int TN, typename Container>
void Matrix<TM, TN, Container>::
initialize(std::initializer_list<real> list) {
	assert_eq(this->SIZE, list.size());
	int i = 0;
	for (auto value : list) {
		this->values[i++] = value;
	}
}


template<int TM, int TN, typename Container>
Matrix<TN, TM, Container> Matrix<TM, TN, Container>::
transpose() const {
	Matrix<TN, TM, Container> result;
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			result(j, i) = (*this)(i, j);
		}
	}
	return result;
}


template<int TM, int TN, typename Container>
void Matrix<TM, TN, Container>::
transposeInplace() {
	(*this) = transpose();
}


template<int TM, int TN, typename Container>
Matrix<TN, TM, Container> Matrix<TM, TN, Container>::
invert() const {
	return GslUtils::invert(*this);
}


template<int TM, int TN, typename Container>
void Matrix<TM, TN, Container>::
invertInplace() {
	(*this) = invert();
}


template<int TM, int TN, typename Container>
Matrix<TM, 1, DefaultMatrixContainer<TM, 1> > Matrix<TM, TN, Container>::
diagonalMultiply(const Matrix<TN, TM>& B) const {
	assert_eq(TM, TN);
	Matrix<TM, 1, DefaultMatrixContainer<TM, 1> > ans;
	for (int i = 0; i < TM; i++) {
		ans(i) = 0;
		for (int j = 0; j < TM; j++) {
			ans(i) += (*this)(i, j) * B(j, i);
		}
	}
	return ans;
}


template<int TM, int TN, typename Container>
real Matrix<TM, TN, Container>::
trace() const {
	assert_eq(TM, TN);
	real ans = 0.0;
	for (int i = 0; i < TN; i++) {
		ans += (*this)(i, i);
	}
	return ans;
}


typedef Matrix<1, 1> Matrix11;
typedef Matrix<2, 2> Matrix22;
typedef Matrix<3, 3> Matrix33;

/**
 * Arbitrary type 2x2 determinant
 */
template<
        typename T11, typename T12,
        typename T21, typename T22
        >
inline auto determinant(const T11 &m11, const T12 &m12,
                        const T21 &m21, const T22 &m22)
->decltype(m11 * m22) {
	return m11 * m22 - m12 * m21;
}

inline real determinant(const Matrix22& m) {
	return determinant(m(0, 0), m(0, 1),
	                   m(1, 0), m(1, 1));
}


/**
 * Arbitrary type 3x3 determinant
 */
template<
        typename T11, typename T12, typename T13,
        typename T21, typename T22, typename T23,
        typename T31, typename T32, typename T33
        >
inline auto determinant(const T11 &m11, const T12 &m12, const T13 &m13,
                        const T21 &m21, const T22 &m22, const T23 &m23,
                        const T31 &m31, const T32 &m32, const T33 &m33)
->decltype(m11 * m22 * m33) {
	return m11 * (m22 * m33 - m23 * m32) -
	       m12 * (m21 * m33 - m23 * m31) +
	       m13 * (m21 * m32 - m22 * m31);
}

inline real determinant(const Matrix33& m) {
	return determinant(m(0, 0), m(0, 1), m(0, 2),
	                   m(1, 0), m(1, 1), m(1, 2),
	                   m(2, 0), m(2, 1), m(2, 2));
}


}
}

namespace std {
template<int TM, int TN, typename Container>
inline std::ostream& operator<<(std::ostream& os,
                                const gcm::linal::Matrix<TM, TN, Container>& matrix) {

	os << std::endl;
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			os << matrix(i, j) << "\t";
		}
		os << "\n";
	}

	return os;
}


}

#endif // LIBGCM_LINAL_MATRIX_HPP
