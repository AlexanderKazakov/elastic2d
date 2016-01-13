#include <gtest/gtest.h>

#include "lib/util/areas/AxisAlignedBoxArea.hpp"
#include "lib/util/areas/SphereArea.hpp"
#include "lib/util/areas/StraightBoundedCylinderArea.hpp"

using namespace gcm;

TEST(Area, AxisAlignedBoxArea) {
	Area* area = new AxisAlignedBoxArea({-5, 0, 2}, {-3, 2, 4});
	ASSERT_EQ(false, area->contains({-10, 1, 3}));
	ASSERT_EQ(false, area->contains({-4, 0, 3}));
	ASSERT_EQ(true, area->contains({-4, 1, 3}));
	ASSERT_ANY_THROW(new AxisAlignedBoxArea({1, 0, 2}, {-3, 2, 4}));
	ASSERT_ANY_THROW(new AxisAlignedBoxArea({-4, 2, 2}, {-3, 2, 4}));
}

TEST(Area, SphereArea) {
	Area* area = new SphereArea(5, {0, 0, 0});
	ASSERT_EQ(false, area->contains({-10, 1, 3}));
	ASSERT_EQ(false, area->contains({5, 0, 0}));
	ASSERT_EQ(true, area->contains({-4, 0, 2.5}));

	area = new SphereArea(1, {1, 2, 3});
	ASSERT_EQ(false, area->contains({0, 1, 5}));
	ASSERT_EQ(true, area->contains({1, 2, 3}));
	ASSERT_EQ(true, area->contains({1.5, 1.5, 2.5}));

	ASSERT_ANY_THROW(new SphereArea(0, {-3, 2, 4}));
	ASSERT_ANY_THROW(new SphereArea(-1.0, {-3, 2, 4}));
}

TEST(Area, StraightBoundedCylinderArea) {
	Area* area = new StraightBoundedCylinderArea(1, {0, 0, 0}, {1, 0, 0});
	ASSERT_EQ(false, area->contains({-10, 1, 3}));
	ASSERT_EQ(false, area->contains({0, 0, 0}));
	ASSERT_EQ(false, area->contains({0.5, 1, 0}));
	ASSERT_EQ(true, area->contains({0.5, 0.5, 0.5}));

	area = new StraightBoundedCylinderArea(2, {1, 1, 1}, {2, 2, 2});
	ASSERT_EQ(false, area->contains({-10, 11, 5}));
	ASSERT_EQ(true, area->contains({1.5, 1.5, 1.5}));

	ASSERT_ANY_THROW(new StraightBoundedCylinderArea(0, {-3, 2, 4}, {-10, 1, 3}));
	ASSERT_ANY_THROW(new StraightBoundedCylinderArea(-1.0, {-3, 2, 4}, {-10, 1, 3}));
	ASSERT_ANY_THROW(new StraightBoundedCylinderArea(1.0, {-3, 2, 4}, {-3, 2, 4}));
}