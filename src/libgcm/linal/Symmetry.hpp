#ifndef LIBGCM_LINAL_SYMMETRY_HPP
#define LIBGCM_LINAL_SYMMETRY_HPP


namespace gcm {
namespace linal {

/** Usual non-symmetric matrix */
struct NonSymmetric { };
/** Symmetric matrix -- quadratic only */
struct Symmetric { };
/** Diagonal matrix -- quadratic only */
struct Diagonal { };


/**
 * A struct for matrix properties connected with its symmetry
 * @tparam TM number of matrix rows
 * @tparam TN number of matrix columns
 */
template<typename TSymmetry, int TM, int TN> struct SymmProps;

template<int TM, int TN>
struct SymmProps<NonSymmetric, TM, TN> {
	static const int M = TM;       ///< number of rows
	static const int N = TN;       ///< number of columns
	static const int SIZE = M * N; ///< size of storage in units of stored elements
	
	static int getIndex(const int i, const int j) {
	/// in such indexation, values in memory placed row by row
		return i * N + j;
	}
};

template<int TM, int TN>
struct SymmProps<Symmetric, TM, TN> {
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

template<int TM, int TN>
struct SymmProps<Diagonal, TM, TN> {
	static const int M = TM;       ///< number of rows
	static const int N = TM;       ///< number of columns
	static const int SIZE = M;     ///< size of storage in units of stored elements
	
	/// This symmetry requires partial MatrixBase specialization, so no getIndex here.
};


/// @name An instrument to specify the symmetry of results 
/// of different arithmetic operations for compiler
/// @{

template<typename TSymmetry1, typename TSymmetry2>
struct LessSymmetric {
	typedef NonSymmetric type;
};

template<>
struct LessSymmetric<Symmetric, Symmetric> {
	typedef Symmetric type;
};

template<>
struct LessSymmetric<Diagonal, Diagonal> {
	typedef Diagonal type;
};

template<>
struct LessSymmetric<Symmetric, Diagonal> {
	typedef Symmetric type;
};

template<>
struct LessSymmetric<Diagonal, Symmetric> {
	typedef Symmetric type;
};

/// @}

} // namespace linal
} // namespace gcm



#endif // LIBGCM_LINAL_SYMMETRY_HPP
