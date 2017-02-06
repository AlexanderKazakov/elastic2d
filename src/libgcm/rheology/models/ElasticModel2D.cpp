#include <libgcm/rheology/models/ElasticModel.hpp>


namespace gcm {


template<>
void ElasticModel<2>::constructNotRotated(
		GcmMatricesPtr m, const real rho,
		const real c11, const real c12, const real c22, const real c66) {
	
	const real cp1 = sqrt(c11/rho);
	const real cp2 = sqrt(c22/rho);
	const real cs  = sqrt(c66/rho);
	
	m->m[0].A = {
		0,0,   -1.0/rho,0,0,
		0,0,    0,-1.0/rho,0,
		
		-c11,0, 0,0,0,
		0,-c66, 0,0,0,
		-c12,0, 0,0,0,
	};
	
	m->m[0].L = {
		-cs, cs, -cp1, cp1, 0,
	};
	
	m->m[0].U = {
			0,   1.0,   0,             1.0/(rho*cs),  0,
			0,   1.0,   0,            -1.0/(rho*cs),  0,
			1.0, 0,     1.0/(rho*cp1), 0,             0,
			1.0, 0,    -1.0/(rho*cp1), 0,             0,
			0,   0,     1.0/c11,         0,          -1.0/c12,
	};
	
	m->m[0].U1 = {
			0,           0,          0.5,            0.5,          0,
			0.5,         0.5,        0,              0,            0,
			0,           0,          0.5*rho*cp1,   -0.5*rho*cp1,  0,
			0.5*rho*cs, -0.5*rho*cs, 0,              0,            0,
			0,           0,          0.5*c12/cp1,   -0.5*c12/cp1, -c12,
	};
	
	
	m->m[1].A = {
		0,0,    0,-1.0/rho,0,
		0,0,    0,0,-1.0/rho,
		
		0,-c12, 0,0,0,
		-c66,0, 0,0,0,
		0,-c22, 0,0,0,
	};
	
	m->m[1].L = {
		-cs, cs, -cp2, cp2, 0,
	};
	
	m->m[1].U = {
			1.0,       0,      0,    1.0/(rho*cs),  0,
			1.0,       0,      0,   -1.0/(rho*cs),  0,
			0,         1.0,    0,    0,             1.0/(rho*cp2),
			0,         1.0,    0,    0,            -1.0/(rho*cp2),
			0,         0,      1.0,  0,            -c12/c22,
	};
	
	m->m[1].U1 = {
			0.5,         0.5,          0,             0,            0,
			0,           0,            0.5,           0.5,          0,
			0,           0,            0.5*c12/cp2,  -0.5*c12/cp2,  1.0,
			0.5*rho*cs, -0.5*rho*cs,   0,             0,            0,
			0,           0,            0.5*rho*cp2,  -0.5*rho*cp2,  0,
	};
	
	m->checkDecomposition();
}


template<>
void ElasticModel<2>::constructGcmMatrices(GcmMatricesPtr m,
		std::shared_ptr<const OrthotropicMaterial> material,
		const MatrixDD& basis) {
	assert_true(basis == linal::identity(basis)); // TODO
	m->basis = basis;
	constructNotRotated(m, material->rho,
			material->c11, material->c12, material->c22, material->c66);
}


template<>
void ElasticModel<1>::constructGcmMatrices(GcmMatricesPtr,
		std::shared_ptr<const OrthotropicMaterial>,
		const MatrixDD&) {
	THROW_UNSUPPORTED("OrthotropicMaterial in 1D is meaningless");
}


}
