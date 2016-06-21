#ifndef LIBGCM_LINAL_SYMMETRY_HPP
#define LIBGCM_LINAL_SYMMETRY_HPP


namespace gcm {
namespace linal {

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


/// @name An instrument to specify the symmetry of results 
/// of different arithmetic operations for compiler
/// @{

template<template<int, int> class TSymmetry1, 
         template<int, int> class TSymmetry2>
struct LessSymmetric {
	template<int TM, int TN>
	using type = NonSymmetric<TM, TN>;
};

template<>
struct LessSymmetric<Symmetric, Symmetric> {
	template<int TM, int TN>
	using type = Symmetric<TM, TN>;
};

template<>
struct LessSymmetric<Diagonal, Diagonal> {
	template<int TM, int TN>
	using type = Diagonal<TM, TN>;
};

template<>
struct LessSymmetric<Symmetric, Diagonal> {
	template<int TM, int TN>
	using type = Symmetric<TM, TN>;
};

template<>
struct LessSymmetric<Diagonal, Symmetric> {
	template<int TM, int TN>
	using type = Symmetric<TM, TN>;
};

/// @}

} // namespace linal
} // namespace gcm



#endif // LIBGCM_LINAL_SYMMETRY_HPP
