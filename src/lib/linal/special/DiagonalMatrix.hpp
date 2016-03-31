#ifndef LIBGCM_LINAL_DIAGONALMATRIX_HPP
#define LIBGCM_LINAL_DIAGONALMATRIX_HPP

namespace gcm {
namespace linal {
/**
 * Special value container for diagonal matrix.
 * @tparam M strings/columns number in matrix
 */
template<int M>
class DiagonalMatrixContainer {
public:
	static const int SIZE = M; ///< size of storage in units of gcm::real
	real values[SIZE];
};


template<int TM>
using DiagonalMatrix = Matrix<TM, TM, DiagonalMatrixContainer<TM> >;

/**
 * Diagonal matrix class (partial specialization)
 * @tparam TM strings/columns number in matrix
 */
template<int TM>
class Matrix<TM, TM, DiagonalMatrixContainer<TM> > : public DiagonalMatrixContainer<TM> {
public:
	static const int M = TM; ///< number of strings
	static const int N = TM; ///< number of columns

	/** Default constructor. */
	Matrix() { }
	/** @param values Values to initialize matrix with, string by string */
	Matrix(std::initializer_list<real> list) {
		this->initialize(list);
	}

	/** @param values Values to initialize matrix with, string by string */
	void initialize(std::initializer_list<real> list);

	/** Read-only access to matrix component */
	real operator()(const int i, const int j) const {
		return (i == j) ? this->values[i] : 0.0;
	}

	/** Read/write access to matrix component */
	real& operator()(const int i, const int j) {
		if (i == j) { return this->values[i]; }
		THROW_INVALID_ARG("Trying to write non-diagonal values of diagonal matrix");
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
	Matrix<TM, TM, DiagonalMatrixContainer<TM> > transpose() const;

	/** Transposes square matrix (modifying matrix contents) */
	void transposeInplace();

	/** @return inverted matrix */
	Matrix<TM, TM, DiagonalMatrixContainer<TM> > invert() const;

	/** Inverses matrix modifying its contents */
	void invertInplace();

	/** @return trace of the matrix */
	real trace() const;

};


template<int TM>
void Matrix<TM, TM, DiagonalMatrixContainer<TM> >::
initialize(std::initializer_list<real> list) {
	assert_eq(this->SIZE, list.size());
	int i = 0;
	for (auto value : list) {
		this->values[i++] = value;
	}
}


template<int TM>
Matrix<TM, TM, DiagonalMatrixContainer<TM> > Matrix<TM, TM, DiagonalMatrixContainer<TM> >::
transpose() const {
	return (*this);
}


template<int TM>
void Matrix<TM, TM, DiagonalMatrixContainer<TM> >::
transposeInplace() { }

template<int TM>
Matrix<TM, TM, DiagonalMatrixContainer<TM> > Matrix<TM, TM, DiagonalMatrixContainer<TM> >::
invert() const {
	Matrix<TM, TM, DiagonalMatrixContainer<TM> > result;
	for (int i = 0; i < TM; i++) {
		assert_ne((*this)(i), 0.0);
		result(i) = 1.0 / (*this)(i);
	}
	return result;
}


template<int TM>
void Matrix<TM, TM, DiagonalMatrixContainer<TM> >::
invertInplace() {
	(*this) = invert();
}


template<int TM>
real Matrix<TM, TM, DiagonalMatrixContainer<TM> >::
trace() const {
	real ans = 0.0;
	for (int i = 0; i < TM; i++) {
		ans += (*this)(i);
	}
	return ans;
}


};
};

#endif // LIBGCM_LINAL_DIAGONALMATRIX_HPP
