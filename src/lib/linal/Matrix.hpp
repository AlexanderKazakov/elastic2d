#ifndef LIBGCM_LINAL_MATRIX_HPP
#define LIBGCM_LINAL_MATRIX_HPP


#include <vector>
#include <iostream>
#include <cmath>

#include <lib/util/Types.hpp>
#include <lib/util/Assertion.hpp>


namespace gcm {
namespace linal {

/**
 * Default implementation for matrix values container.
 * @tparam TSize size of storage
 * @tparam TElement type of stored elements
 */
template<int TSize, typename TElement>
struct DefaultContainer {

	typedef TElement ElementType;
	
	/// @name Construction and assignment @{
	
	DefaultContainer() = default;
	
	DefaultContainer(const DefaultContainer& orig) = default;
	
	DefaultContainer(DefaultContainer&& orig) = default;
	
	DefaultContainer(const std::initializer_list<TElement>& list) {
		assert_eq(TSize, list.size());
		std::copy(list.begin(), list.end(), values);
	}
	
	DefaultContainer(const std::vector<TElement>& list) {
		assert_eq(TSize, list.size());
		std::copy(list.begin(), list.end(), values);
	}
	
	template<typename TElement2>
	DefaultContainer(const DefaultContainer<TSize, TElement2>& orig) {
		for (int i = 0; i < TSize; i++) {
			(*this)(i) = static_cast<TElement>(orig(i));
		}
	}
	
	template<typename TElement2>
	DefaultContainer(DefaultContainer<TSize, TElement2>&& orig) {
		for (int i = 0; i < TSize; i++) {
			(*this)(i) = static_cast<TElement>(orig(i));
		}
	}
	

	DefaultContainer& operator=(const DefaultContainer& m2) = default;
	
	DefaultContainer& operator=(DefaultContainer&& m2) = default;
	
	template<typename TElement2>
	DefaultContainer& operator=(const DefaultContainer<TSize, TElement2>& orig) {
		return (*this) = DefaultContainer(orig);
	}

	template<typename TElement2>
	DefaultContainer& operator=(DefaultContainer<TSize, TElement2>&& orig) {
		return (*this) = DefaultContainer(orig);
	}
	
	void copyFrom(const DefaultContainer& origin) {
		(*this) = origin;
	}
	
	/// @}
	
	/// @name Access to array elements @{
	/** Read-only access */
	TElement operator()(const int i) const {
		return values[i];
	}

	/** Read/write access */
	TElement& operator()(const int i) {
		return values[i];
	}
	///@}
		
private:
	TElement values[TSize]; ///< the storage
};


/**
 * Usual non-symmetric matrix properties
 * @tparam TM number of rows
 * @tparam TN number of columns
 */
template<int TM, int TN>
struct NonSymmetric {
	static const int M = TM;       ///< number of rows
	static const int N = TN;       ///< number of columns
	static const int SIZE = M * N; ///< size of storage in units of stored elements
	
	static int getIndex(const int i, const int j) {
	/// in such indexation, values in memory placed row by row
		return i * N + j;
	}
};


/**
 * Symmetric matrix properties. TM == TN.
 * @tparam TM number of rows
 * @tparam TN number of columns
 */
template<int TM, int TN = TM>
struct Symmetric {
	static const int M = TM;                 ///< number of rows
	static const int N = TM;                 ///< number of columns
	static const int SIZE = M * (M + 1) / 2; ///< size of storage in units of stored elements
	
	static int getIndex(const int i, const int j) {
	/// in such indexation, values in memory placed row by row,
	/// upper part only (first row - M elements, second - (M-1), etc)
		static_assert(TM == TN, "Symmetric matrix is square");
		return (i < j) ? i * N - ((i - 1) * i) / 2 + j - i
		               : j * N - ((j - 1) * j) / 2 + i - j;
	}
};


/**
 * Diagonal matrix properties. TM == TN.
 * @tparam TM number of rows
 * @tparam TN number of columns
 */
template<int TM, int TN = TM>
struct Diagonal {
	static const int M = TM;       ///< number of rows
	static const int N = TM;       ///< number of columns
	static const int SIZE = M;     ///< size of storage in units of stored elements
	
	/// This symmetry requires partial MatrixBase specialization, so no getIndex here.
};


/**
 * Base class for matrices.
 * @tparam TM number of rows
 * @tparam TN number of columns
 * @tparam TElement type of stored elements
 * @tparam TSymmetry NonSymmetric, Symmetric, Diagonal, etc ..
 * @tparam TContainer type of storage
 */
template<int TM, int TN,
         typename TElement,
         template<int, int> class TSymmetry,
         template<int, typename> class TContainer>
struct MatrixBase :
		TContainer<TSymmetry<TM, TN>::SIZE, TElement> {
		
	typedef TSymmetry<TM, TN> Symmetry;
	static const int  M   = Symmetry::M;
	static const int  N   = Symmetry::N;
	static const int SIZE = Symmetry::SIZE;

	typedef MatrixBase<M, 1, TElement, NonSymmetric, TContainer>  ColumnType;
	typedef MatrixBase<1, N, TElement, NonSymmetric, TContainer>  RowType;
	
	typedef TContainer<SIZE, TElement>                            Container;

	using Container::Container;
	using Container::operator();
	
	/** @return zero-filled matrix */
	inline static MatrixBase Zeros();
	/** @return matrix with all ones */
	inline static MatrixBase Ones();
	/** @return identity matrix */
	inline static MatrixBase Identity();

	///@name Access to matrix elements @{
	/** Read-only access */
	TElement operator()(const int i, const int j) const {
		return (*this)(Symmetry::getIndex(i, j));
	}

	/** Read/write access */
	TElement& operator()(const int i, const int j) {
		return (*this)(Symmetry::getIndex(i, j));
	}
	///@}

	///@name Access to rows and columns @{
	/** @return j-th column. */
	ColumnType getColumn(const int j) const {
		ColumnType ans;
		for (int i = 0; i < M; i++) {
			ans(i) = (*this)(i, j);
		}
		return ans;
	}

	/** set j-th column */
	void setColumn(const int j, const ColumnType& column) {
		for (int i = 0; i < M; i++) {
			(*this)(i, j) = column(i);
		}
	}
	
	/** @return i-th row. */
	RowType getRow(const int i) const {
		RowType ans;
		for (int j = 0; j < N; j++) {
			ans(j) = (*this)(i, j);
		}
		return ans;
	}
	
	/** set i-th row */
	void setRow(const int i, const RowType& row) {
		for (int j = 0; j < N; j++) {
			(*this)(i, j) = row(j);
		}
	}
	///@}
	
	///@name Operations inplace (rewrite the object with operation result) @{
	void transposeInplace() { (*this) = transpose(*this); }
	void invertInplace()    { (*this) = invert(*this); }
	///@}
};


/**
 * Base for diagonal matrices.
 * @tparam TM number of rows == number of columns
 * @tparam TElement type of stored elements
 * @tparam TContainer type of storage
 */
template<int TM,
         typename TElement,
         template<int, typename> class TContainer>
struct MatrixBase<TM, TM, TElement, Diagonal, TContainer> :
		TContainer<Diagonal<TM, TM>::SIZE, TElement> {
		
	typedef Diagonal<TM, TM> Symmetry;
	static const int  M   = Symmetry::M;
	static const int  N   = Symmetry::N;
	static const int SIZE = Symmetry::SIZE;

	typedef MatrixBase<M, 1, TElement, NonSymmetric, TContainer>  ColumnType;
	typedef MatrixBase<1, N, TElement, NonSymmetric, TContainer>  RowType;
	
	typedef TContainer<SIZE, TElement>                            Container;

	using Container::Container;
	using Container::operator();

	/** @return zero-filled matrix */
	inline static MatrixBase Zeros();
	/** @return matrix with all ones */
	inline static MatrixBase Ones();
	/** @return identity matrix */
	inline static MatrixBase Identity();

	///@name Access to elements @{
	/** 
	 * Read-only access to matrix components.
	 * For read/write access use operator()(const int) from container.
	 */
	TElement operator()(const int i, const int j) const {
		return (i == j) ? (*this)(i) : 0;
	}
	///@}
};


/// @name Different special matrices derived from MatrixBase
/// @{

template<int TM, int TN> struct Matrix :
	MatrixBase<TM, TN, real, NonSymmetric, DefaultContainer> {
	
	typedef MatrixBase<TM, TN, real, NonSymmetric, DefaultContainer> Base;
	using Base::Base;
};
	
template<int TM> struct Vector :
	Matrix<TM, 1> {
	
	typedef Matrix<TM, 1> Base;
	using Base::Base;
};

template<int TN> struct Row :
	Matrix<1, TN> {
	
	typedef Matrix<1, TN> Base;
	using Base::Base;
};


template<int TM> struct SymmetricMatrix :
	MatrixBase<TM, TM, real, Symmetric, DefaultContainer> {
	
	typedef MatrixBase<TM, TM, real, Symmetric, DefaultContainer> Base;
	using Base::Base;
};


template<int TM> struct DiagonalMatrix :
	MatrixBase<TM, TM, real, Diagonal, DefaultContainer> {
	
	typedef MatrixBase<TM, TM, real, Diagonal, DefaultContainer> Base;
	using Base::Base;
};


template<int TM, int TN> struct MatrixInt :
	MatrixBase<TM, TN, int, NonSymmetric, DefaultContainer> {
	
	typedef MatrixBase<TM, TN, int, NonSymmetric, DefaultContainer> Base;
	using Base::Base;
};

template<int TM> struct VectorInt :
	MatrixInt<TM, 1> {
	
	typedef MatrixInt<TM, 1> Base;
	using Base::Base;
};


typedef Matrix<1, 1> Matrix11;
typedef Matrix<2, 2> Matrix22;
typedef Matrix<3, 3> Matrix33;

typedef Vector<1> Real1;
typedef Vector<2> Real2;
typedef Vector<3> Real3;
typedef Vector<4> Real4;

typedef VectorInt<1> Int1;
typedef VectorInt<2> Int2;
typedef VectorInt<3> Int3;


template<int TM, int TN, typename TElement> struct MATRIX :
	MatrixBase<TM, TN, TElement, NonSymmetric, DefaultContainer> {
	
	typedef MatrixBase<TM, TN, TElement, NonSymmetric, DefaultContainer> Base;
	using Base::Base;
};
	
template<int TM, typename TElement> struct VECTOR :
	MATRIX<TM, 1, TElement> {
	
	typedef MATRIX<TM, 1, TElement> Base;
	using Base::Base;
};

template<int TM, typename TElement> struct SYMMETRIC_MATRIX :
	MatrixBase<TM, TM, TElement, Symmetric, DefaultContainer> {
	
	typedef MatrixBase<TM, TM, TElement, Symmetric, DefaultContainer> Base;
	using Base::Base;
};

/// @}


/// @name Fill in the matrix with some values
/// @{

/**
 * Bottom of recursion for zeros function for matrices
 */
template<typename TNumber>
typename std::enable_if<std::is_arithmetic<TNumber>::value, TNumber>::type
zeros(const TNumber&) {
	return 0; 
}

/**
 * @return zero-filled matrix
 */
template<int TM, int TN,
         typename TElement,
         template<int, int> class TSymmetry,
         template<int, typename> class TContainer>
MatrixBase<TM, TN, TElement, TSymmetry, TContainer>
zeros(const MatrixBase<TM, TN, TElement, TSymmetry, TContainer>&) {
	MatrixBase<TM, TN, TElement, TSymmetry, TContainer> result;
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			result(i, j) = zeros(result(i, j));
		}
	}
	return result;
}

/**
 * @return zero-filled diagonal matrix
 */
template<int TM,
         typename TElement,
         template<int, typename> class TContainer>
MatrixBase<TM, TM, TElement, Diagonal, TContainer>
zeros(const MatrixBase<TM, TM, TElement, Diagonal, TContainer>&) {
	MatrixBase<TM, TM, TElement, Diagonal, TContainer> result;
	for (int i = 0; i < TM; i++) {
		result(i) = zeros(result(i));
	}
	return result;
}


/**
 * Bottom of recursion for ones function for matrices
 */
template<typename TNumber>
typename std::enable_if<std::is_arithmetic<TNumber>::value, TNumber>::type
ones(const TNumber&) {
	return 1; 
}

/**
 * @return matrix with all unit elements
 */
template<int TM, int TN,
         typename TElement,
         template<int, int> class TSymmetry,
         template<int, typename> class TContainer>
MatrixBase<TM, TN, TElement, TSymmetry, TContainer>
ones(const MatrixBase<TM, TN, TElement, TSymmetry, TContainer>&) {
	MatrixBase<TM, TN, TElement, TSymmetry, TContainer> result;
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			result(i, j) = ones(result(i, j));
		}
	}
	return result;
}

/**
 * @return diagonal matrix with all unit elements at diagonal
 */
template<int TM,
         typename TElement,
         template<int, typename> class TContainer>
MatrixBase<TM, TM, TElement, Diagonal, TContainer>
ones(const MatrixBase<TM, TM, TElement, Diagonal, TContainer>&) {
	MatrixBase<TM, TM, TElement, Diagonal, TContainer> result;
	for (int i = 0; i < TM; i++) {
		result(i) = ones(result(i));
	}
	return result;
}


/**
 * Bottom of recursion for identity function for matrices
 */
template<typename TNumber>
typename std::enable_if<std::is_arithmetic<TNumber>::value, TNumber>::type
identity(const TNumber&) {
	return 1; 
}

/**
 * @return identity (aka unity) matrix
 */
template<int TM,
         typename TElement,
         template<int, int> class TSymmetry,
         template<int, typename> class TContainer>
MatrixBase<TM, TM, TElement, TSymmetry, TContainer>
identity(const MatrixBase<TM, TM, TElement, TSymmetry, TContainer>&) {
	MatrixBase<TM, TM, TElement, TSymmetry, TContainer> result;
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TM; j++) {
			if (i == j) {
				result(i, j) = identity(result(i, j));		
			} else {
				result(i, j) = zeros(result(i, j));
			}
		}
	}
	return result;
}

/**
 * @return identity (aka unity) matrix with TSymmetry == Diagonal
 */
template<int TM,
         typename TElement,
         template<int, typename> class TContainer>
MatrixBase<TM, TM, TElement, Diagonal, TContainer>
identity(const MatrixBase<TM, TM, TElement, Diagonal, TContainer>& m) {
	return ones(m);
}


/**
 * Fill in the matrix with zeros
 */
template<int TM, int TN,
         typename TElement,
         template<int, int> class TSymmetry,
         template<int, typename> class TContainer>
MatrixBase<TM, TN, TElement, TSymmetry, TContainer>&
clear(MatrixBase<TM, TN, TElement, TSymmetry, TContainer>& m) {
	return (m = zeros(m));
}
///@}


template<int TM, int TN,
         typename TElement,
         template<int, int> class TSymmetry,
         template<int, typename> class TContainer>
inline
MatrixBase<TM, TN, TElement, TSymmetry, TContainer>
MatrixBase<TM, TN, TElement, TSymmetry, TContainer>::
Zeros() {
	return zeros(MatrixBase<TM, TN, TElement, TSymmetry, TContainer>());
}

template<int TM,
         typename TElement,
         template<int, typename> class TContainer>
inline
MatrixBase<TM, TM, TElement, Diagonal, TContainer>
MatrixBase<TM, TM, TElement, Diagonal, TContainer>::
Zeros() {
	return zeros(MatrixBase<TM, TM, TElement, Diagonal, TContainer>());
}


template<int TM, int TN,
         typename TElement,
         template<int, int> class TSymmetry,
         template<int, typename> class TContainer>
inline
MatrixBase<TM, TN, TElement, TSymmetry, TContainer>
MatrixBase<TM, TN, TElement, TSymmetry, TContainer>::
Ones() {
	return ones(MatrixBase<TM, TN, TElement, TSymmetry, TContainer>());
}

template<int TM,
         typename TElement,
         template<int, typename> class TContainer>
inline
MatrixBase<TM, TM, TElement, Diagonal, TContainer>
MatrixBase<TM, TM, TElement, Diagonal, TContainer>::
Ones() {
	return ones(MatrixBase<TM, TM, TElement, Diagonal, TContainer>());
}


template<int TM, int TN,
         typename TElement,
         template<int, int> class TSymmetry,
         template<int, typename> class TContainer>
inline
MatrixBase<TM, TN, TElement, TSymmetry, TContainer>
MatrixBase<TM, TN, TElement, TSymmetry, TContainer>::
Identity() {
	static_assert(TM == TN, "");
	return identity(MatrixBase<TM, TM, TElement, TSymmetry, TContainer>());
}

template<int TM,
         typename TElement,
         template<int, typename> class TContainer>
inline
MatrixBase<TM, TM, TElement, Diagonal, TContainer>
MatrixBase<TM, TM, TElement, Diagonal, TContainer>::
Identity() {
	return identity(MatrixBase<TM, TM, TElement, Diagonal, TContainer>());
}


}
}



namespace std {

template<int TM, int TN,
         typename TElement,
         template<int, int> class TSymmetry,
         template<int, typename> class TContainer>
inline std::ostream& operator<<(std::ostream& os, 
		const gcm::linal::MatrixBase<TM, TN, TElement, TSymmetry, TContainer>& matrix) {

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
