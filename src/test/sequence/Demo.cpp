#include <gtest/gtest.h>

TEST(Demo, sizes) {
	struct Empty { };

	struct A {
		double a[2];
		Empty empty;
	};

	struct B : public Empty {
		double b[2];
	};

	std::cout << "sizeof(A) = " << sizeof(A) << std::endl; // 24
	std::cout << "sizeof(B) = " << sizeof(B) << std::endl; // 16

	std::cout << "sizeof(A[1000]) = " << sizeof(A[1000]) << std::endl; // 24000
	std::cout << "sizeof(B[1000]) = " << sizeof(B[1000]) << std::endl; // 16000
}