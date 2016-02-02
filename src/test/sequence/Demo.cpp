#include <gtest/gtest.h>
#include <memory>

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

	std::cout << "sizeof(A[1000000]) = " << sizeof(A[1000000]) << std::endl; // 24000000
	std::cout << "sizeof(B[1000000]) = " << sizeof(B[1000000]) << std::endl; // 16000000

	std::cout << "sizeof(std::shared_ptr<A>) = " << sizeof(std::shared_ptr<A>) << std::endl; // 16
	std::cout << "sizeof(A*) = " << sizeof(A*) << std::endl; // 8
}