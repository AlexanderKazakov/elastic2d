#include <gtest/gtest.h>

#include <libgcm/util/math/Area.hpp>

using namespace gcm;

TEST(Area, InfiniteArea) {
	Area* area = new InfiniteArea();
	ASSERT_EQ(true, area->contains({-10, 1, 300}));
	ASSERT_EQ(true, area->contains({-4e+55, 2e+11, 3e-99}));
	ASSERT_EQ(true, area->contains({0, 0, 0}));

	area->move({1e+200, 1, 2});
	ASSERT_EQ(true, area->contains({-10, 1, 3}));
	ASSERT_EQ(true, area->contains({-4e+55, 2e+11, 3e-99}));
	ASSERT_EQ(true, area->contains({0, 0, 0}));
}

TEST(Area, AxisAlignedBoxArea) {
	Area* area = new AxisAlignedBoxArea({-5, 0, 2}, {-3, 2, 4});
	ASSERT_EQ(false, area->contains({-10, 1, 3}));
	ASSERT_EQ(false, area->contains({-4, 0, 3}));
	ASSERT_EQ(true, area->contains({-4, 1, 3}));

	area->move({1, 1, 1});
	ASSERT_EQ(false, area->contains({-9, 2, 4}));
	ASSERT_EQ(false, area->contains({-3, 1, 4}));
	ASSERT_EQ(true, area->contains({-3, 2, 4}));

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

	area->move({5, 5, -5});
	ASSERT_EQ(false, area->contains({5, 6, 0}));
	ASSERT_EQ(true, area->contains({6, 7, -2}));
	ASSERT_EQ(true, area->contains({6.5, 6.5, -2.5}));

	ASSERT_ANY_THROW(new SphereArea(0, {-3, 2, 4}));
	ASSERT_ANY_THROW(new SphereArea(-1.0, {-3, 2, 4}));
}

TEST(Area, StraightBoundedCylinderArea) {
	Area* area = new StraightBoundedCylinderArea(1, {0, 0, 0}, {1, 0, 0});
	ASSERT_EQ(false, area->contains({-10, 1, 3}));
	ASSERT_EQ(false, area->contains({0, 0, 0}));
	ASSERT_EQ(false, area->contains({0.5, 1, 0}));
	ASSERT_EQ(true, area->contains({0.5, 0.5, 0.5}));

	area->move({-10, 0, 10});
	ASSERT_EQ(false, area->contains({-20, 1, 13}));
	ASSERT_EQ(false, area->contains({-10, 0, 10}));
	ASSERT_EQ(false, area->contains({-9.5, 1, 10}));
	ASSERT_EQ(true, area->contains({-9.5, 0.5, 10.5}));

	area = new StraightBoundedCylinderArea(2, {1, 1, 1}, {2, 2, 2});
	ASSERT_EQ(false, area->contains({-10, 11, 5}));
	ASSERT_EQ(true, area->contains({1.5, 1.5, 1.5}));

	ASSERT_ANY_THROW(new StraightBoundedCylinderArea(0, {-3, 2, 4}, {-10, 1, 3}));
	ASSERT_ANY_THROW(new StraightBoundedCylinderArea(-1.0, {-3, 2, 4}, {-10, 1, 3}));
	ASSERT_ANY_THROW(new StraightBoundedCylinderArea(1.0, {-3, 2, 4}, {-3, 2, 4}));
}
