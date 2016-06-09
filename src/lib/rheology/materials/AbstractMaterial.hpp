#ifndef LIBGCM_ABSTRACTMATERIAL_HPP
#define LIBGCM_ABSTRACTMATERIAL_HPP

#include <initializer_list>

#include <lib/util/Enum.hpp>
#include <lib/linal/linal.hpp>

namespace gcm {
struct AbstractMaterial {

	/** 6x6 symmetric elastic matrix */
	typedef linal::SymmetricMatrix<6> ElasticMatrix;
	
	/** 
	 * 3^4 elastic tensor.
	 * Symmetry q(i, j)(k, l) = q(j, i)(k, l) = q(i, j)(l, k) is already
	 * provided by symmetry of matrices.
	 * Symmetry q(i, j)(k, l) = q(k, l)(i, j) must be provided by user.
	 * @see convert
	 */
	typedef linal::SYMMETRIC_MATRIX<3, linal::SymmetricMatrix<3>> ElasticTensor;
	
	virtual ~AbstractMaterial() = default;

	/** Convert 3^4 elastic tensor to 6x6 elastic matrix */
	static ElasticMatrix convert(const ElasticTensor& q) {
		auto ans = ElasticMatrix::Zeros();
	
		ans(0, 0) = q(0, 0)(0, 0); 
		ans(1, 1) = q(1, 1)(1, 1); 
		ans(2, 2) = q(2, 2)(2, 2); 
		ans(0, 1) = q(0, 0)(1, 1); 
		ans(0, 2) = q(0, 0)(2, 2); 
		ans(1, 2) = q(1, 1)(2, 2); 
		ans(0, 3) = q(0, 0)(1, 2); 
		ans(0, 4) = q(0, 0)(0, 2); 
		ans(0, 5) = q(0, 0)(0, 1); 
		ans(1, 3) = q(1, 1)(1, 2); 
		ans(1, 4) = q(1, 1)(0, 2); 
		ans(1, 5) = q(1, 1)(0, 1); 
		ans(2, 3) = q(2, 2)(1, 2); 
		ans(2, 4) = q(2, 2)(0, 2); 
		ans(2, 5) = q(2, 2)(0, 1); 
		ans(4, 4) = q(0, 2)(0, 2); 
		ans(3, 3) = q(1, 2)(1, 2); 
		ans(5, 5) = q(0, 1)(0, 1); 
		ans(3, 4) = q(0, 2)(1, 2); 
		ans(3, 5) = q(0, 1)(1, 2); 
		ans(4, 5) = q(0, 1)(0, 2);
		
		return ans;
	}
	
	
	/** Convert 6x6 elastic matrix to 3^4 elastic tensor */
	static ElasticTensor convert(const ElasticMatrix& c) {
		auto q = ElasticTensor::Zeros();
		
		q(0, 0)(0, 0) = c(0, 0); 
		q(1, 1)(1, 1) = c(1, 1); 
		q(2, 2)(2, 2) = c(2, 2); 
	
		q(1, 2)(1, 2) = c(3, 3); 
		q(0, 2)(0, 2) = c(4, 4); 
		q(0, 1)(0, 1) = c(5, 5); 
		
		q(0, 0)(1, 1) = q(1, 1)(0, 0) = c(0, 1); 
		q(0, 0)(2, 2) = q(2, 2)(0, 0) = c(0, 2); 
		q(1, 1)(2, 2) = q(2, 2)(1, 1) = c(1, 2); 
	
		q(0, 0)(1, 2) = q(1, 2)(0, 0) = c(0, 3); 
		q(0, 0)(0, 2) = q(0, 2)(0, 0) = c(0, 4); 
		q(0, 0)(0, 1) = q(0, 1)(0, 0) = c(0, 5); 
	
		q(1, 1)(1, 2) = q(1, 2)(1, 1) = c(1, 3); 
		q(1, 1)(0, 2) = q(0, 2)(1, 1) = c(1, 4); 
		q(1, 1)(0, 1) = q(0, 1)(1, 1) = c(1, 5); 
	
		q(2, 2)(1, 2) = q(1, 2)(2, 2) = c(2, 3); 
		q(2, 2)(0, 2) = q(0, 2)(2, 2) = c(2, 4); 
		q(2, 2)(0, 1) = q(0, 1)(2, 2) = c(2, 5); 
	
		q(1, 2)(0, 2) = q(0, 2)(1, 2) = c(3, 4); 
		q(1, 2)(0, 1) = q(0, 1)(1, 2) = c(3, 5); 
		q(0, 2)(0, 1) = q(0, 1)(0, 2) = c(4, 5);
		
		return q;
	}
	
	
	/**
	 * Return given tensor in system of axes which is rotated 
	 * firstly, by phi(0) radians around x-axis clockwise,
	 * secondly, by phi(1) radians around y-axis clockwise,
	 * thirdly, by phi(2) radians around z-axis clockwise
	 */
	static ElasticTensor rotate(const ElasticTensor& t, const Real3& phi) {
		auto ans = ElasticTensor::Zeros();
		
		linal::Matrix33 G = linal::getZRotationMatrix(phi(2)) *
		                    linal::getYRotationMatrix(phi(1)) * 
		                    linal::getXRotationMatrix(phi(0));
		
		for (int m = 0; m < 3; m++)
		for (int n = m; n < 3; n++) 
		for (int p = 0; p < 3; p++) 
		for (int q = p; q < 3; q++) {
			
			for (int i = 0; i < 3; i++) 
			for (int j = 0; j < 3; j++) 
			for (int k = 0; k < 3; k++) 
			for (int l = 0; l < 3; l++) {
				ans(m, n)(p, q) += 
						G(m, i) * G(n, j) * G(p, k) * G(q, l) * t(i, j)(k, l);
			}
		}
		
		return ans;
	}


	/**
	 * Return given elastic matrix in system of axes which is rotated 
	 * firstly, by phi(0) radians around x-axis clockwise,
	 * secondly, by phi(1) radians around y-axis clockwise,
	 * thirdly, by phi(2) radians around z-axis clockwise
	 */
	static ElasticMatrix rotate(const ElasticMatrix& c, const Real3& phi) {
		return convert(rotate(convert(c), phi));
	}
	
	
};


}

#endif // LIBGCM_ABSTRACTMATERIAL_HPP
