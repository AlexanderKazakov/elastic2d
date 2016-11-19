#include <gtest/gtest.h>
#include <memory>
#include <cmath>

#include <libgcm/linal/linal.hpp>

//#include <libgcm/rheology/materials/OrthotropicMaterial.hpp>
using namespace gcm;

//TEST(Demo, sizes) {
//	struct Empty { };

//	struct A {
//		double a[2];
//		Empty empty;
//	};

//	struct B : public Empty {
//		double b[2];
//	};

//	struct Cb {
//		double c[2];
//		virtual void foo() = 0;

//	};

//	struct C : public Cb {
//		virtual void foo() override { }
//	};

//	struct Db {
//		double d[2];
//		virtual void foo() = 0;
//		virtual void foo2() = 0;
//		virtual void foo3() = 0;
//		virtual void foo4() = 0;
//		virtual void foo5() = 0;
//		virtual int foo6() = 0;

//	};

//	struct D : public Db {
//		virtual void foo() override { }
//		virtual void foo2() override { }
//		virtual void foo3() override { }
//		virtual void foo4() override { }
//		virtual void foo5() override { }
//		virtual int foo6() override { return 0; }
//	};
	
///*	
//		std::cout << sizeof(long double) << std::endl;
//		std::cout << sizeof(double) << std::endl;

//		std::cout << "sizeof(A) = " << sizeof(A) << std::endl; // 24
//	   std::cout << "sizeof(B) = " << sizeof(B) << std::endl; // 16

//	   std::cout << "sizeof(A[1000000]) = " << sizeof(A[1000000]) << std::endl; // 24000000
//	   std::cout << "sizeof(B[1000000]) = " << sizeof(B[1000000]) << std::endl; // 16000000

//	   std::cout << "sizeof(std::shared_ptr<A>) = " << sizeof(std::shared_ptr<A>) << std::endl;
//	      // 16
//	   std::cout << "sizeof(A*) = " << sizeof(A*) << std::endl; // 8

//	   std::cout << "sizeof(Cb) = " << sizeof(Cb) << std::endl; // 24
//	   std::cout << "sizeof(C) = " << sizeof(C) << std::endl; // 24
//	   std::cout << "sizeof(C[1000000]) = " << sizeof(C[1000000]) << std::endl; // 24000000

//	   std::cout << "sizeof(Db) = " << sizeof(Db) << std::endl; // 24
//	   std::cout << "sizeof(D) = " << sizeof(D) << std::endl; // 24
//	   std::cout << "sizeof(D[1000000]) = " << sizeof(D[1000000]) << std::endl; // 24000000
//	 */
//}


//TEST(Random, rand) {
//	Utils::seedRand();
//	for (int i = 0; i < 10; i++) {
//		std::cout << linal::randomBasis(linal::Matrix22()) << std::endl;
//	}
//}

/*TEST(OrthotropicMaterial, rotate) {
	auto material = std::make_shared<OrthotropicMaterial>(
			OrthotropicMaterial(1, {11, 12, 13, 22, 23, 33, 44, 55, 66}));
	
	auto m = material->getElasticMatrix();
	std::cout.precision(2);
	std::cout << m << std::endl
			<< OrthotropicMaterial::rotate(m, {M_PI / 4, M_PI / 3, M_PI / 6}) << std::endl;
}*/

