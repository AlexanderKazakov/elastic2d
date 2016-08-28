#include <gtest/gtest.h>

#include <libgcm/util/StringUtils.hpp>

using namespace gcm;


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
