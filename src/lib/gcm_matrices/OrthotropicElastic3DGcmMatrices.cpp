#include "lib/gcm_matrices/OrthotropicElastic3DGcmMatrices.hpp"

using namespace gcm;

const std::map<Waves::T, int/* number of column in U1 */> OrthotropicElastic3DGcmMatrices::WAVE_COLUMNS = {
		{Waves::T::P_FORWARD,   5},
		{Waves::T::P_BACKWARD,  4},
		{Waves::T::S1_FORWARD,  1},
		{Waves::T::S1_BACKWARD, 0},
		{Waves::T::S2_FORWARD,  3},
		{Waves::T::S2_BACKWARD, 2}
};

OrthotropicElastic3DGcmMatrices::OrthotropicElastic3DGcmMatrices(const real &rho, const real &lambda, const real &mu) :
		rho(rho) /* TODO Material */ {

	// TODO Material
	c11 =  90000;
	c22 = 180000;
	c33 = 360000;

	c12 =  70000;
	c13 =  50000;
	c23 =  40000;

	c44 =  60000;
	c55 =  30000;
	c66 =  20000;

	m[0].A = linal::Matrix<9, 9>({0, 0, 0, -1.0 / rho, 0, 0, 0, 0, 0,
	                              0, 0, 0, 0, -1.0 / rho, 0, 0, 0, 0,
	                              0, 0, 0, 0, 0, -1.0 / rho, 0, 0, 0,
	                              -c11, 0, 0, 0, 0, 0, 0, 0, 0,
	                              0, -c66, 0, 0, 0, 0, 0, 0, 0,
	                              0, 0, -c55, 0, 0, 0, 0, 0, 0,
	                              -c12, 0, 0, 0, 0, 0, 0, 0, 0,
	                              0, 0, 0, 0, 0, 0, 0, 0, 0,
	                              -c13, 0, 0, 0, 0, 0, 0, 0, 0});

	m[0].L.createDiagonal(
			{-sqrt(c66 / rho), sqrt(c66 / rho), -sqrt(c55 / rho), sqrt(c55 / rho), -sqrt(c11 / rho), sqrt(c11 / rho), 0, 0, 0});

	m[0].U = linal::Matrix<9, 9>({0, 1.0, 0, 0, 1.0 / (sqrt(c66) * sqrt(rho)), 0, 0, 0, 0,
	                              0, 1.0, 0, 0, -1.0 / (sqrt(c66) * sqrt(rho)), 0, 0, 0, 0,
	                              0, 0, 1.0, 0, 0, 1.0 / (sqrt(c55) * sqrt(rho)), 0, 0, 0,
	                              0, 0, 1.0, 0, 0, -1.0 / (sqrt(c55) * sqrt(rho)), 0, 0, 0,
	                              1.0, 0, 0, 1.0 / (sqrt(c11) * sqrt(rho)), 0, 0, 0, 0, 0,
	                              1.0, 0, 0, -1.0 / (sqrt(c11) * sqrt(rho)), 0, 0, 0, 0, 0,
	                              0, 0, 0, -c12 / c11, 0, 0, 1.0, 0, 0.0,
	                              0, 0, 0, 0, 0, 0, 0, 1.0, 0,
	                              0, 0, 0, -(1.0 * c13) / c11, 0, 0, 0, 0, 1.0});

	m[0].U1 = linal::Matrix<9, 9>({0, 0, 0, 0, 0.5, 0.5, 0, 0, 0,
	                               0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
	                               0, 0, 0.5, 0.5, 0, 0, 0, 0, 0,
	                               0, 0, 0, 0, 0.5 * sqrt(c11) * sqrt(rho), -0.5 * sqrt(c11) * sqrt(rho), 0, 0, 0,
	                               0.5 * sqrt(c66) * sqrt(rho), -0.5 * sqrt(c66) * sqrt(rho), 0, 0, 0, 0, 0, 0, 0,
	                               0, 0, 0.5 * sqrt(c55) * sqrt(rho), -0.5 * sqrt(c55) * sqrt(rho), 0, 0, 0, 0, 0,
	                               0, 0, 0, 0, (0.5 * c12 * sqrt(rho)) / sqrt(c11),
	                               -(0.5 * c12 * sqrt(rho)) / sqrt(c11), 1, 0, 0,
	                               0, 0, 0, 0, 0, 0, 0, 1, 0,
	                               0, 0, 0, 0, (0.5 * c13 * sqrt(rho)) / sqrt(c11), -(0.5 * c13 * sqrt(rho)) / sqrt(c11), 0, 0, 1});


	m[1].A = linal::Matrix<9, 9>({0, 0, 0, 0, -1.0 / rho, 0, 0, 0, 0,
	                              0, 0, 0, 0, 0, 0, -1.0 / rho, 0, 0,
	                              0, 0, 0, 0, 0, 0, 0, -1.0 / rho, 0,
	                              0, -c12, 0, 0, 0, 0, 0, 0, 0,
	                              -c66, 0, 0, 0, 0, 0, 0, 0, 0,
	                              0, 0, 0, 0, 0, 0, 0, 0, 0,
	                              0, -c22, 0, 0, 0, 0, 0, 0, 0,
	                              0, 0, -c44, 0, 0, 0, 0, 0, 0,
	                              0, -c23, 0, 0, 0, 0, 0, 0, 0});

	m[1].L.createDiagonal(
			{-sqrt(c66 / rho), sqrt(c66 / rho), -sqrt(c44 / rho), sqrt(c44 / rho), -sqrt(c22 / rho), sqrt(c22 / rho), 0, 0, 0});

	m[1].U = linal::Matrix<9, 9>({1.0, 0, 0, 0, 1.0 / (sqrt(c66) * sqrt(rho)), 0, 0, 0, 0,
	                              1.0, 0, 0, 0, -1.0 / (sqrt(c66) * sqrt(rho)), 0, 0, 0, 0,
	                              0, 0, 1.0, 0, 0, 0, 0, 1.0 / (sqrt(c44) * sqrt(rho)), 0,
	                              0, 0, 1.0, 0, 0, 0, 0, -1.0 / (sqrt(c44) * sqrt(rho)), 0,
	                              0, 1.0, 0, 0, 0, 0, 1.0 / (sqrt(c22) * sqrt(rho)), 0, 0,
	                              0, 1.0, 0, 0, 0, 0, -1.0 / (sqrt(c22) * sqrt(rho)), 0, 0,
	                              0, 0, 0, 1.0, 0, 0, -c12 / c22, 0, 0.0,
	                              0, 0, 0, 0, 0, 1.0, 0, 0, 0,
	                              0, 0, 0, 0, 0, 0, -(1.0 * c23) / c22, 0, 1.0});

	m[1].U1 = linal::Matrix<9, 9>({0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
	                               0, 0, 0, 0, 0.5, 0.5, 0, 0, 0,
	                               0, 0, 0.5, 0.5, 0, 0, 0, 0, 0,
	                               0, 0, 0, 0, (0.5 * c12) / sqrt(c22 / rho), -(0.5 * c12) / sqrt(c22 / rho), 1, 0, 0,
	                               0.5 * sqrt(c66) * sqrt(rho), -0.5 * sqrt(c66) * sqrt(rho), 0, 0, 0, 0, 0, 0, 0,
	                               0, 0, 0, 0, 0, 0, 0, 1, 0,
	                               0, 0, 0, 0, 0.5 * sqrt(c22) * sqrt(rho), -0.5 * sqrt(c22) * sqrt(rho), 0, 0, 0,
	                               0, 0, 0.5 * sqrt(c44) * sqrt(rho), -0.5 * sqrt(c44) * sqrt(rho), 0, 0, 0, 0, 0,
	                               0, 0, 0, 0, (0.5 * c23 * sqrt(rho)) / sqrt(c22), -(0.5 * c23 * sqrt(rho)) / sqrt(c22), 0, 0, 1});


	m[2].A = linal::Matrix<9, 9>({0, 0, 0, 0, 0, -1.0 / rho, 0, 0, 0,
	                              0, 0, 0, 0, 0, 0, 0, -1.0 / rho, 0,
	                              0, 0, 0, 0, 0, 0, 0, 0, -1.0 / rho,
	                              0, 0, -c13, 0, 0, 0, 0, 0, 0,
	                              0, 0, 0, 0, 0, 0, 0, 0, 0,
	                              -c55, 0, 0, 0, 0, 0, 0, 0, 0,
	                              0, 0, -c23, 0, 0, 0, 0, 0, 0,
	                              0, -c44, 0, 0, 0, 0, 0, 0, 0,
	                              0, 0, -c33, 0, 0, 0, 0, 0, 0});

	m[2].L.createDiagonal(
			{-sqrt(c55 / rho), sqrt(c55 / rho), -sqrt(c44 / rho), sqrt(c44 / rho), -sqrt(c33 / rho), sqrt(c33 / rho), 0, 0, 0});

	m[2].U = linal::Matrix<9, 9>({1.0, 0, 0, 0, 0, 1.0 / (sqrt(c55) * sqrt(rho)), 0, 0, 0,
	                              1.0, 0, 0, 0, 0, -1.0 / (sqrt(c55) * sqrt(rho)), 0, 0, 0,
	                              0, 1.0, 0, 0, 0, 0, 0, 1.0 / (sqrt(c44) * sqrt(rho)), 0,
	                              0, 1.0, 0, 0, 0, 0, 0, -1.0 / (sqrt(c44) * sqrt(rho)), 0,
	                              0, 0, 1.0, 0, 0, 0, 0, 0, 1.0 / (sqrt(c33) * sqrt(rho)),
	                              0, 0, 1.0, 0, 0, 0, 0, 0, -1.0 / (sqrt(c33) * sqrt(rho)),
	                              0, 0, 0, 1.0, 0, 0, 0, 0, -(1.0 * c13) / c33,
	                              0, 0, 0, 0, 1.0, 0, 0, 0, 0,
	                              0, 0, 0, 0, 0, 0, 1.0, 0, -(1.0 * c23) / c33});

	m[2].U1 = linal::Matrix<9, 9>({0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
	                               0, 0, 0.5, 0.5, 0, 0, 0, 0, 0,
	                               0, 0, 0, 0, 0.5, 0.5, 0, 0, 0,
	                               0, 0, 0, 0, (0.5 * c13 * sqrt(rho)) / sqrt(c33),
	                               -(0.5 * c13 * sqrt(rho)) / sqrt(c33), 1, 0, 0,
	                               0, 0, 0, 0, 0, 0, 0, 1, 0,
	                               0.5 * sqrt(c55) * sqrt(rho), -0.5 * sqrt(c55) * sqrt(rho), 0, 0, 0, 0, 0, 0, 0,
	                               0, 0, 0, 0, (0.5 * c23 * sqrt(rho)) / sqrt(c33), -(0.5 * c23 * sqrt(rho)) / sqrt(c33), 0, 0, 1,
	                               0, 0, 0.5 * sqrt(c44) * sqrt(rho), -0.5 * sqrt(c44) * sqrt(rho), 0, 0, 0, 0, 0,
	                               0, 0, 0, 0, 0.5 * sqrt(c33) * sqrt(rho), -0.5 * sqrt(c33) * sqrt(rho), 0, 0, 0});

}