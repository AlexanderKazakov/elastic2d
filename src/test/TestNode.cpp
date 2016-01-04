#include <gtest/gtest.h>

#include "lib/config.hpp"
#include "lib/nodes/IdealElastic2DNode.hpp"

using namespace gcm;

TEST(IdealElastic2DNode, Access)
{
	IdealElastic2DNode node;
	node(0) = 0; node(1) = 1; node(2) = 2; node(3) = 3; node(4) = 4;
	ASSERT_EQ(node.Vx, 0);
	ASSERT_EQ(node.Vy, 1);
	ASSERT_EQ(node.Sxx, 2);
	ASSERT_EQ(node.Sxy, 3);
	ASSERT_EQ(node.Syy, 4);

	IdealElastic2DNode node2;
	node2.Vx = 0; node2.Vy = 1;
	node2.Sxx = 2; node2.Sxy = 3; node2.Syy = 4;
	ASSERT_EQ(node(0), 0);
	ASSERT_EQ(node(1), 1);
	ASSERT_EQ(node(2), 2);
	ASSERT_EQ(node(3), 3);
	ASSERT_EQ(node(4), 4);

}
