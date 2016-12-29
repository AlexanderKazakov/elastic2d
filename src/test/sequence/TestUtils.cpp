#include <gtest/gtest.h>

#include <libgcm/util/StringUtils.hpp>
#include <libgcm/util/math/Histogram.hpp>

using namespace gcm;


TEST(Histogram, hist) {
	std::vector<real> a(1001);
	for (size_t i = 0; i < a.size(); i++) { a[i] = (real)i; }
	Histogram ha(a.begin(), a.end(), a.size());
	ASSERT_EQ(0, ha.min());
	ASSERT_EQ(a.size() - 1, ha.max());
	ASSERT_NEAR(real(a.size() - 1) / 2, ha.mean(), EQUALITY_TOLERANCE);
	ASSERT_EQ(std::vector<size_t>(a.size(), 1), ha.binCounts());
	real binSize = real(a.size() - 1) / real(a.size());
	std::vector<real> centers = ha.binCenters();
	for (size_t i = 0; i < centers.size(); i++) {
		ASSERT_NEAR((real(i) + 0.5) * binSize, centers[i], EQUALITY_TOLERANCE);
	}
}


TEST(StringUtils, split) {
	const std::string s = "1223 32sjd op   34";
	const auto v = StringUtils::split(s, ' ');
	ASSERT_EQ(4, v.size());
	ASSERT_EQ("1223",  v[0]);
	ASSERT_EQ("32sjd", v[1]);
	ASSERT_EQ("op",    v[2]);
	ASSERT_EQ("34",    v[3]);
	
	const auto v2 = StringUtils::split(s, '2');
	ASSERT_EQ(3, v2.size());
	ASSERT_EQ("1",     v2[0]);
	ASSERT_EQ("3 3",   v2[1]);
	ASSERT_EQ("sjd op   34", v2[2]);
	
	const auto v3 = StringUtils::split(s, '\t');
	ASSERT_EQ(1, v3.size());
	ASSERT_EQ("1223 32sjd op   34", v3[0]);
}


TEST(StringUtils, toString) {
	ASSERT_EQ("0000000004", StringUtils::toString(4, 10));
	ASSERT_EQ("1234", StringUtils::toString(1234, 4));
}
