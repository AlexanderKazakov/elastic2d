#include <cmath>

#include "GcmMatrices2D.hpp"

using namespace gcm;
using namespace gcm::linal;

GcmMatrices2D::GcmMatrices2D(const real& rho, const real& lambda, const real& mu) :
		rho(rho), lambda(lambda), mu(mu) {

	Ax.A(0, 0) = 0;              Ax.A(0, 1) = 0;              Ax.A(0, 2) = -1.0/rho; Ax.A(0, 3) = 0;        Ax.A(0, 4) = 0;
	Ax.A(1, 0) = 0;              Ax.A(1, 1) = 0;              Ax.A(1, 2) = 0;        Ax.A(1, 3) = -1.0/rho; Ax.A(1, 4) = 0;
	Ax.A(2, 0) = -(lambda+2*mu); Ax.A(2, 1) = 0;              Ax.A(2, 2) = 0;        Ax.A(2, 3) = 0;        Ax.A(2, 4) = 0;
	Ax.A(3, 0) = 0;              Ax.A(3, 1) = -mu;            Ax.A(3, 2) = 0;        Ax.A(3, 3) = 0;        Ax.A(3, 4) = 0;
	Ax.A(4, 0) = -lambda;        Ax.A(4, 1) = 0;              Ax.A(4, 2) = 0;        Ax.A(4, 3) = 0;        Ax.A(4, 4) = 0;

	Ay.A(0, 0) = 0;              Ay.A(0, 1) = 0;              Ay.A(0, 2) = 0;        Ay.A(0, 3) = -1.0/rho; Ay.A(0, 4) = 0;
	Ay.A(1, 0) = 0;              Ay.A(1, 1) = 0;              Ay.A(1, 2) = 0;        Ay.A(1, 3) = 0;        Ay.A(1, 4) = -1.0/rho;
	Ay.A(2, 0) = 0;              Ay.A(2, 1) = -lambda;        Ay.A(2, 2) = 0;        Ay.A(2, 3) = 0;        Ay.A(2, 4) = 0;
	Ay.A(3, 0) = -mu;            Ay.A(3, 1) = 0;              Ay.A(3, 2) = 0;        Ay.A(3, 3) = 0;        Ay.A(3, 4) = 0;
	Ay.A(4, 0) = 0;              Ay.A(4, 1) = -(lambda+2*mu); Ay.A(4, 2) = 0;        Ay.A(4, 3) = 0;        Ay.A(4, 4) = 0;


	Ax.L.createDiagonal( {-sqrt((lambda+2*mu)/rho),sqrt((lambda+2*mu)/rho),-sqrt(mu/rho),sqrt(mu/rho),0} );

	Ax.U = Matrix( {1.0,0,1.0/(sqrt(rho)*sqrt(lambda+2*mu)),0,0,
	                1.0,0,-1.0/(sqrt(rho)*sqrt(lambda+2*mu)),0,0,
	                0,1.0,0,1.0/(sqrt(mu)*sqrt(rho)),0,
	                0,1.0,0,-1.0/(sqrt(mu)*sqrt(rho)),0,
	                0,0,1.0/(lambda+2*mu),0,-1.0/lambda} );

	Ax.U1 = Matrix( {0.5,0.5,0,0,0,
	                 0,0,0.5,0.5,0,
	                 0.5*sqrt(rho)*sqrt(lambda+2*mu),-0.5*sqrt(rho)*sqrt(lambda+2*mu),0,0,0,
	                 0,0,0.5*sqrt(mu)*sqrt(rho),-0.5*sqrt(mu)*sqrt(rho),0,
	                 (0.5*sqrt(rho)*lambda)/sqrt(lambda+2*mu),-(0.5*sqrt(rho)*lambda)/sqrt(lambda+2*mu),0,0,-lambda} );

	Ay.L.createDiagonal( {-sqrt((lambda+2*mu)/rho),sqrt((lambda+2*mu)/rho),-sqrt(mu/rho),sqrt(mu/rho),0} );

	Ay.U = Matrix( {0,1.0,0,0,1.0/(sqrt(rho)*sqrt(lambda+2*mu)),
	                0,1.0,0,0,-1.0/(sqrt(rho)*sqrt(lambda+2*mu)),
	                1.0,0,0,1.0/(sqrt(mu)*sqrt(rho)),0,
	                1.0,0,0,-1.0/(sqrt(mu)*sqrt(rho)),0,
	                0,0,1.0,0,-(1.0*lambda)/(lambda+2*mu)} );

	Ay.U1 = Matrix( {0,0,0.5,0.5,0,
	                 0.5,0.5,0,0,0,
	                 (0.5*lambda)/sqrt((lambda+2*mu)/rho),-(0.5*lambda)/sqrt((lambda+2*mu)/rho),0,0,1.0,
	                 0,0,0.5*sqrt(mu)*sqrt(rho),-0.5*sqrt(mu)*sqrt(rho),0,
	                 0.5*sqrt(rho)*sqrt(lambda+2*mu),-0.5*sqrt(rho)*sqrt(lambda+2*mu),0,0,0} );
}


const GcmMatrix &GcmMatrices2D::A(const int stage) const {
	if      (stage == 0) return Ax;
	else if (stage == 1) return Ay;
	else throw "Invalid stage number";
}