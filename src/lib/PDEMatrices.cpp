#include <string.h>
#include <cmath>

#include "lib/PDEMatrices.hpp"


PDEMatrices::PDEMatrices(const real& rho, const real& lambda, const real& mu) :
		rho(rho), lambda(lambda), mu(mu) {

	bool variablesOrderIsCorrect = (
			static_cast<uint>(NodeMap::Vx) == 0 &&
			static_cast<uint>(NodeMap::Vy) == 1 &&
			static_cast<uint>(NodeMap::Sxx) == 2 &&
			static_cast<uint>(NodeMap::Sxy) == 3 &&
			static_cast<uint>(NodeMap::Syy) == 4 );

	if ( !variablesOrderIsCorrect ) throw "Variables order is not correct";

	Ax.A(0, 0) = 0;              Ax.A(0, 1) = 0;   Ax.A(0, 2) = -1.0/rho; Ax.A(0, 3) = 0;        Ax.A(0, 4) = 0;
	Ax.A(1, 0) = 0;              Ax.A(1, 1) = 0;   Ax.A(1, 2) = 0;        Ax.A(1, 3) = -1.0/rho; Ax.A(1, 4) = 0;
	Ax.A(2, 0) = -(lambda+2*mu); Ax.A(2, 1) = 0;   Ax.A(2, 2) = 0;        Ax.A(2, 3) = 0;        Ax.A(2, 4) = 0;
	Ax.A(3, 0) = 0;              Ax.A(3, 1) = -mu; Ax.A(3, 2) = 0;        Ax.A(3, 3) = 0;        Ax.A(3, 4) = 0;
	Ax.A(4, 0) = -lambda;        Ax.A(4, 1) = 0;   Ax.A(4, 2) = 0;        Ax.A(4, 3) = 0;        Ax.A(4, 4) = 0;

	Ay.A(0, 0) = 0;   Ay.A(0, 1) = 0;              Ay.A(0, 2) = 0; Ay.A(0, 3) = -1.0/rho; Ay.A(0, 4) = 0;
	Ay.A(1, 0) = 0;   Ay.A(1, 1) = 0;              Ay.A(1, 2) = 0; Ay.A(1, 3) = 0;        Ay.A(1, 4) = -1.0/rho;
	Ay.A(2, 0) = 0;   Ay.A(2, 1) = -lambda;        Ay.A(2, 2) = 0; Ay.A(2, 3) = 0;        Ay.A(2, 4) = 0;
	Ay.A(3, 0) = -mu; Ay.A(3, 1) = 0;              Ay.A(3, 2) = 0; Ay.A(3, 3) = 0;        Ay.A(3, 4) = 0;
	Ay.A(4, 0) = 0;   Ay.A(4, 1) = -(lambda+2*mu); Ay.A(4, 2) = 0; Ay.A(4, 3) = 0;        Ay.A(4, 4) = 0;


	Ax.L.createDiagonal({-sqrt((lambda+2*mu)/rho),sqrt((lambda+2*mu)/rho),-sqrt(mu/rho),sqrt(mu/rho),0});

	Ax.U.createMatrix( { {1.0,0,1.0/(sqrt(rho)*sqrt(lambda+2*mu)),0,0},
	                     {1.0,0,-1.0/(sqrt(rho)*sqrt(lambda+2*mu)),0,0},
	                     {0,1.0,0,1.0/(sqrt(mu)*sqrt(rho)),0},
	                     {0,1.0,0,-1.0/(sqrt(mu)*sqrt(rho)),0},
	                     {0,0,1.0/(lambda+2*mu),0,-1.0/lambda} } );

	Ax.U1.createMatrix( { {0.5,0.5,0,0,0},
	                      {0,0,0.5,0.5,0},
	                      {0.5*sqrt(rho)*sqrt(lambda+2*mu),-0.5*sqrt(rho)*sqrt(lambda+2*mu),0,0,0},
	                      {0,0,0.5*sqrt(mu)*sqrt(rho),-0.5*sqrt(mu)*sqrt(rho),0},
	                      {(0.5*sqrt(rho)*lambda)/sqrt(lambda+2*mu),-(0.5*sqrt(rho)*lambda)/sqrt(lambda+2*mu),0,0,-lambda} } );

	Ay.L.createDiagonal({-sqrt((lambda+2*mu)/rho),sqrt((lambda+2*mu)/rho),-sqrt(mu/rho),sqrt(mu/rho),0});

	Ay.U.createMatrix( { {0,1.0,0,0,1.0/(sqrt(rho)*sqrt(lambda+2*mu))},
	                     {0,1.0,0,0,-1.0/(sqrt(rho)*sqrt(lambda+2*mu))},
	                     {1.0,0,0,1.0/(sqrt(mu)*sqrt(rho)),0},
	                     {1.0,0,0,-1.0/(sqrt(mu)*sqrt(rho)),0},
	                     {0,0,1.0,0,-(1.0*lambda)/(lambda+2*mu)} } );

	Ay.U1.createMatrix( { {0,0,0.5,0.5,0},
	                      {0.5,0.5,0,0,0},
	                      {(0.5*lambda)/sqrt((lambda+2*mu)/rho),-(0.5*lambda)/sqrt((lambda+2*mu)/rho),0,0,1.0},
	                      {0,0,0.5*sqrt(mu)*sqrt(rho),-0.5*sqrt(mu)*sqrt(rho),0},
	                      {0.5*sqrt(rho)*sqrt(lambda+2*mu),-0.5*sqrt(rho)*sqrt(lambda+2*mu),0,0,0} } );
}


const PDEMatrix &PDEMatrices::A(const uint stage) const {
	if      (stage == 0) return Ax;
	else if (stage == 1) return Ay;
	else throw "Invalid stage number";
}


void Vector::createVector(const std::initializer_list<real> &list) {
	if(list.size() != N) throw "Size of initializer_list isn't equal to N";

	uint i = 0;
	for (auto& value : list) {
		(*this)(i) = value;
		i++;
	}
}


Vector Vector::operator*(const real &b) const {
	Vector ans;
	for (uint i = 0; i < N; i++) {
		ans(i) = get(i) * b;
	}
	return ans;
}


Vector Vector::operator-(const Vector &b) const {
	Vector ans;
	for (uint i = 0; i < N; i++) {
		ans(i) = get(i) - b.get(i);
	}
	return ans;
}


void Vector::operator+=(const Vector &b) {
	for (uint i = 0; i < N; i++) {
		(*this)(i) += b.get(i);
	}
}


bool Vector::operator==(const Vector &b) const {
	for (uint i = 0; i < N; i++) {
		if (fabs(get(i) - b.get(i)) > EQUALITY_TOLERANCE) return false;
	}
	return true;
}


void Matrix::createMatrix(const std::initializer_list<std::initializer_list<real>> &list) {
	if(list.size() != N) throw "Size of initializer_list isn't equal to N";

	uint i = 0;
	for(auto& str : list) {
		if(str.size() != N) throw "Size of initializer_list isn't equal to N";

		uint j = 0;
		for (auto &value : str) {
			(*this)(i, j) = value;
			j++;
		}
		i++;
	}
}


void Matrix::createDiagonal(const std::initializer_list<real>& list) {
	memset(m, 0, N * N * sizeof(real));
	if(list.size() != N) throw "Size of initializer_list isn't equal to N";

	uint i = 0;
	for (auto& value : list) {
		(*this)(i, i) = value;
		i++;
	}
}


void Matrix::setColumn(const uint i, const Vector &column) {
	for (uint j = 0; j < N; j++) {
		(*this)(j, i) = column.get(j);
	}
}


Vector Matrix::getColumn(const uint i) const {
	Vector ans;
	for (uint j = 0; j < N; j++) {
		ans(j) = get(j, i);
	}
	return ans;
}


Vector Matrix::getDiagonalMultipliedBy(const real &c) const {
	Vector ans;
	for (uint i = 0; i < N; i++) {
		ans(i) = get(i, i) * c;
	}
	return ans;
}


Vector Matrix::diagonalMultiply(const Matrix &B) const {
	Vector ans;
	for (uint i = 0; i < N; i++) {
		ans(i) = 0;
		for (uint j = 0; j < N; j++) {
			ans(i) += get(i, j) * B.get(j, i);
		}
	}
	return ans;
}


Matrix Matrix::operator*(const Matrix &B) const {
	Matrix C; // C = this * B
	for (uint i = 0; i < N; i++) {
		for (uint j = 0; j < N; j++) {
			C(i, j) = 0.0;
			for (uint k = 0; k < N; k++) {
				C(i, j) += get(i, k) * B.get(k, j);
			}
		}
	}
	return C;
}


Vector Matrix::operator*(const Vector &b) const {
	Vector c; // c = this * b
	for (uint i = 0; i < N; i++) {
		c(i) = 0.0;
		for (uint j = 0; j < N; j++) {
			c(i) += get(i, j) * b.get(j);
		}
	}
	return c;
}


real Matrix::getTrace() const {
	real tr = 0.0;
	for (uint i = 0; i < N; i++) {
		tr += get(i, i);
	}
	return tr;
}


