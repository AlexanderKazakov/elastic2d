#include "IdealElastic2DGcmMatrices.hpp"

using namespace gcm;

const std::map<Waves::WAVE, int /* number of column in U1 */> IdealElastic2DGcmMatrices::WAVE_COLUMNS = {
		{Waves::WAVE::P_FORWARD,   1},
		{Waves::WAVE::P_BACKWARD,  0},
		{Waves::WAVE::S1_FORWARD,  3},
		{Waves::WAVE::S1_BACKWARD, 2}
};

IdealElastic2DGcmMatrices::IdealElastic2DGcmMatrices(const real& rho, const real& lambda, const real& mu) :
		rho(rho), lambda(lambda), mu(mu) {

	m[0].A(0, 0) = 0;              m[0].A(0, 1) = 0;              m[0].A(0, 2) = -1.0/rho; m[0].A(0, 3) = 0;        m[0].A(0, 4) = 0;
	m[0].A(1, 0) = 0;              m[0].A(1, 1) = 0;              m[0].A(1, 2) = 0;        m[0].A(1, 3) = -1.0/rho; m[0].A(1, 4) = 0;
	m[0].A(2, 0) = -(lambda+2*mu); m[0].A(2, 1) = 0;              m[0].A(2, 2) = 0;        m[0].A(2, 3) = 0;        m[0].A(2, 4) = 0;
	m[0].A(3, 0) = 0;              m[0].A(3, 1) = -mu;            m[0].A(3, 2) = 0;        m[0].A(3, 3) = 0;        m[0].A(3, 4) = 0;
	m[0].A(4, 0) = -lambda;        m[0].A(4, 1) = 0;              m[0].A(4, 2) = 0;        m[0].A(4, 3) = 0;        m[0].A(4, 4) = 0;

	m[1].A(0, 0) = 0;              m[1].A(0, 1) = 0;              m[1].A(0, 2) = 0;        m[1].A(0, 3) = -1.0/rho; m[1].A(0, 4) = 0;
	m[1].A(1, 0) = 0;              m[1].A(1, 1) = 0;              m[1].A(1, 2) = 0;        m[1].A(1, 3) = 0;        m[1].A(1, 4) = -1.0/rho;
	m[1].A(2, 0) = 0;              m[1].A(2, 1) = -lambda;        m[1].A(2, 2) = 0;        m[1].A(2, 3) = 0;        m[1].A(2, 4) = 0;
	m[1].A(3, 0) = -mu;            m[1].A(3, 1) = 0;              m[1].A(3, 2) = 0;        m[1].A(3, 3) = 0;        m[1].A(3, 4) = 0;
	m[1].A(4, 0) = 0;              m[1].A(4, 1) = -(lambda+2*mu); m[1].A(4, 2) = 0;        m[1].A(4, 3) = 0;        m[1].A(4, 4) = 0;


	m[0].L.createDiagonal( {-sqrt((lambda+2*mu)/rho),sqrt((lambda+2*mu)/rho),-sqrt(mu/rho),sqrt(mu/rho),0} );

	m[0].U = linal::Matrix<5,5>( {1.0,0,1.0/(sqrt(rho)*sqrt(lambda+2*mu)),0,0,
	                              1.0,0,-1.0/(sqrt(rho)*sqrt(lambda+2*mu)),0,0,
	                              0,1.0,0,1.0/(sqrt(mu)*sqrt(rho)),0,
	                              0,1.0,0,-1.0/(sqrt(mu)*sqrt(rho)),0,
	                              0,0,1.0/(lambda+2*mu),0,-1.0/lambda} );

	m[0].U1 = linal::Matrix<5,5>( {0.5,0.5,0,0,0,
	                               0,0,0.5,0.5,0,
	                               0.5*sqrt(rho)*sqrt(lambda+2*mu),-0.5*sqrt(rho)*sqrt(lambda+2*mu),0,0,0,
	                               0,0,0.5*sqrt(mu)*sqrt(rho),-0.5*sqrt(mu)*sqrt(rho),0,
	                               (0.5*sqrt(rho)*lambda)/sqrt(lambda+2*mu),-(0.5*sqrt(rho)*lambda)/sqrt(lambda+2*mu),0,0,-lambda} );

	m[1].L.createDiagonal( {-sqrt((lambda+2*mu)/rho),sqrt((lambda+2*mu)/rho),-sqrt(mu/rho),sqrt(mu/rho),0} );

	m[1].U = linal::Matrix<5,5>( {0,1.0,0,0,1.0/(sqrt(rho)*sqrt(lambda+2*mu)),
	                              0,1.0,0,0,-1.0/(sqrt(rho)*sqrt(lambda+2*mu)),
	                              1.0,0,0,1.0/(sqrt(mu)*sqrt(rho)),0,
	                              1.0,0,0,-1.0/(sqrt(mu)*sqrt(rho)),0,
	                              0,0,1.0,0,-(1.0*lambda)/(lambda+2*mu)} );

	m[1].U1 = linal::Matrix<5,5>( {0,0,0.5,0.5,0,
	                               0.5,0.5,0,0,0,
	                               (0.5*lambda)/sqrt((lambda+2*mu)/rho),-(0.5*lambda)/sqrt((lambda+2*mu)/rho),0,0,1.0,
	                               0,0,0.5*sqrt(mu)*sqrt(rho),-0.5*sqrt(mu)*sqrt(rho),0,
	                               0.5*sqrt(rho)*sqrt(lambda+2*mu),-0.5*sqrt(rho)*sqrt(lambda+2*mu),0,0,0} );
}