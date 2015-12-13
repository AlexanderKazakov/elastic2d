#ifndef ELASTIC2D_MATRIX_HPP
#define ELASTIC2D_MATRIX_HPP

#include <iostream>

#include "lib/config.hpp"

class Vector {
	real v[N];
public:
	inline real& operator()(const int i) {
		return v[i];
	};
	inline const real& get(const int i) const {
		return v[i];
	};
	/**
	 * @param list of values
	 */
	void createVector(const std::initializer_list<real>& list);
	Vector operator*(const real& b) const;
	Vector operator-(const Vector& b) const;
	void operator+=(const Vector& b);
	bool operator==(const Vector& b) const;
};

class Matrix {
	real m[N * N];
public:
	inline real& operator()(const int i, const int j) {
		return m[i * N + j];
	};
	inline const real& get(const int i, const int j) const {
		return m[i * N + j];
	};
	/**
	 * @param list list of lists, where each one represents a matrix string
	 */
	void createMatrix(const std::initializer_list<std::initializer_list<real>>& list);
	/**
	 * @param list of diagonal values
	 */
	void createDiagonal(const std::initializer_list<real>& list);
	/**
	 * Fill in the i-th column with %columns' values
	 */
	void setColumn(const int i, const Vector& column);
	/**
	 * @return i-th column
	 */
	Vector getColumn(const int i) const;
	/**
	 * @return Vector with values from matrix diagonal,
	 * multiplied by %c
	 */
	Vector getDiagonalMultipliedBy(const real& c) const;
	/**
	 * @return Vector composed of diagonal elements of multiplication this matrix by given matrix
	 */
	Vector diagonalMultiply(const Matrix& B) const;
	Matrix operator*(const Matrix& B) const;
	Vector operator*(const Vector& b) const;
	real getTrace() const;
};

/**
 * Represents a matrix in PDE along some direction with its eigensystem
 */
class PDEMatrix {
public:
	Matrix A; // matrix for some axis in PDE
	/* A = U1 * L * U  and so right eigenvectors are columns of the U1 */
	Matrix U; // matrix of left eigenvectors (aka eigenstrings)
	Matrix U1; // matrix of right eigenvectors
	Matrix L; // diagonal eigenvalue matrix
};

/**
 * The PDE is:
 * d\vec{u}/dt + \mat{Ax} \partial{\vec{u}}/\partial{x} + \mat{Ay} \partial{\vec{u}}/\partial{y} = \vec{f}.
 * The class represents \mat{Ax} and \mat{Ay} with their eigensystems for a concrete material.
 */
class PDEMatrices {

	PDEMatrix Ax;
	PDEMatrix Ay;

public:

	real rho;
	real lambda;
	real mu;

	PDEMatrices(const real& rho, const real& lambda, const real& mu);
	/**
	 * @return PDEMatrix along stage direction
	 */
	const PDEMatrix& A(const int stage) const;

};

namespace std {
	inline ostream& operator<<(std::ostream &os, const Vector& vector) {

		os << "Vector:\n";
		for (int i = 0; i < N; i++) {
				os << vector.get(i) << "\n";
		}

		return os;
	};

	inline ostream& operator<<(std::ostream &os, const Matrix& matrix) {

		os << "Matrix:\n";
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				os << matrix.get(i, j) << "\t";
			}
			os << "\n";
		}

		return os;
	};
}

#endif //ELASTIC2D_MATRIX_HPP
