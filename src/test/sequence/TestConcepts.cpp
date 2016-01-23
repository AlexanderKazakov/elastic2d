#include <gtest/gtest.h>

#include <lib/util/Concepts.hpp>

using namespace gcm;

template<class TConcept>
class TestConcepts : public testing::Test {
protected:
	void testSizeOfMap() {
		ASSERT_EQ(static_cast<int>(TConcept::T::SIZE_OF_ENUM), TConcept::NAME.size());
	};
};

/** Look at https://github.com/google/googletest/blob/master/googletest/samples/sample6_unittest.cc for explaination */
#if GTEST_HAS_TYPED_TEST_P
using testing::Types;
TYPED_TEST_CASE_P(TestConcepts);
TYPED_TEST_P(TestConcepts, SizeOfMap) {
	this->testSizeOfMap();
}
REGISTER_TYPED_TEST_CASE_P(TestConcepts, SizeOfMap);

// When new concept arrears, place it here
typedef Types<PhysicalQuantities, Waves, BorderCondition> AllImplementations;

INSTANTIATE_TYPED_TEST_CASE_P(AllConcepts, TestConcepts, AllImplementations);
#endif // GTEST_HAS_TYPED_TEST_P

