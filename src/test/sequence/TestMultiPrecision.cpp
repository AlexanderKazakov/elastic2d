#include <gtest/gtest.h>
#include <libgcm/util/math/MultiPrecision.hpp>

using namespace gcm;


TEST(MultiPrecision, miscellaneous) {
	static_assert(std::is_arithmetic<MultiPrecision>::value, 
			"MultiPrecision must be cpp-arithmetic type, fix namespace std");
	
	std::string sqrt2 = 
			std::string("1.41421356237309504880168872420969807856967187537") +
			std::string("6948073176679737990732478462107038850387534327641");

	MultiPrecision calc = sqrt(MultiPrecision::create("2", sqrt2.length()));
	std::string calcStr = calc.toString().substr(0, sqrt2.size());
	ASSERT_EQ(sqrt2, calcStr);
	
	
	const size_t precision = 50;
	MultiPrecision::setDefaultDecimalPrecision(precision);
	MultiPrecision::Matrix<3,3> A = {
			3, 0, 0,
			0, 5, 0,
			0, 0, 6,
	};
	MultiPrecision::Vector<3> b = {1, 1, 1};
	auto x = linal::solveLinearSystem(A, b);
	ASSERT_EQ("0.333333333333333333333333333333333333333333333333", 
			x(0).toString().substr(0, precision));
	ASSERT_EQ("0.2", 
			x(1).setPrecision(precision).toString().substr(0, precision));
	ASSERT_EQ("0.166666666666666666666666666666666666666666666666", 
			x(2).toString().substr(0, precision));
	
	
	MultiPrecision m1 = MultiPrecision(1) / 3;
	MultiPrecision m4 = MultiPrecision(1) / 70;
	MultiPrecision::Matrix<2, 2> m = {m1, 0, 0, m4};
	MultiPrecision::Vector<2> v = {3, 7};
	MultiPrecision::Vector<2> mv = m * v;
	ASSERT_EQ("1", mv(0).toString().substr(0, precision));
	ASSERT_EQ("0.1", mv(1).setDecimalPrecision(precision - 2).toString().substr(0, precision));
}
