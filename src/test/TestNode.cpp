#include <gtest/gtest.h>

#include "lib/config.hpp"
#include "lib/Node.hpp"


TEST(Node, Access)
{
	Node node;
	node.u(0) = 0; node.u(1) = 1; node.u(2) = 2; node.u(3) = 3; node.u(4) = 4;
	ASSERT_EQ(node.get(NodeMap::Vx), 0);
	ASSERT_EQ(node.get(NodeMap::Vy), 1);
	ASSERT_EQ(node.get(NodeMap::Sxx), 2);
	ASSERT_EQ(node.get(NodeMap::Sxy), 3);
	ASSERT_EQ(node.get(NodeMap::Syy), 4);

	Node node2;
	node2(NodeMap::Vx) = 0; node2(NodeMap::Vy) = 1;
	node2(NodeMap::Sxx) = 2; node2(NodeMap::Sxy) = 3; node2(NodeMap::Syy) = 4;
	ASSERT_EQ(node.u.get(0), 0);
	ASSERT_EQ(node.u.get(1), 1);
	ASSERT_EQ(node.u.get(2), 2);
	ASSERT_EQ(node.u.get(3), 3);
	ASSERT_EQ(node.u.get(4), 4);

}
