#include <gtest/gtest.h>

#include "lib/nodes/IdealElastic3DNode.hpp"
#include "lib/nodes/IdealElastic2DNode.hpp"
#include "lib/nodes/IdealElastic1DNode.hpp"

using namespace gcm;

TEST(IdealElastic3DNode, Access)
{
	IdealElastic3DNode node;
	for (int i = 0; i < node.M; i++) {
		node(i) = i;
	}
	ASSERT_EQ(node.Vx, 0);
	ASSERT_EQ(node.Vy, 1);
	ASSERT_EQ(node.Vz, 2);
	ASSERT_EQ(node.Sxx, 3);
	ASSERT_EQ(node.Sxy, 4);
	ASSERT_EQ(node.Sxz, 5);
	ASSERT_EQ(node.Syy, 6);
	ASSERT_EQ(node.Syz, 7);
	ASSERT_EQ(node.Szz, 8);

	IdealElastic3DNode node2;
	node2.Vx = 0; node2.Vy = 1; node2.Vz = 2;
	node2.Sxx = 3; node2.Sxy = 4; node2.Sxz = 5; node2.Syy = 6; node2.Syz = 7; node2.Szz = 8;
	for (int j = 0; j < node2.M; j++) {
		ASSERT_EQ(node(j), j);
	}
}

TEST(IdealElastic2DNode, Access)
{
	IdealElastic2DNode node;
	for (int i = 0; i < node.M; i++) {
		node(i) = i;
	}
	ASSERT_EQ(node.Vx, 0);
	ASSERT_EQ(node.Vy, 1);
	ASSERT_EQ(node.Sxx, 2);
	ASSERT_EQ(node.Sxy, 3);
	ASSERT_EQ(node.Syy, 4);

	IdealElastic2DNode node2;
	node2.Vx = 0; node2.Vy = 1;
	node2.Sxx = 2; node2.Sxy = 3; node2.Syy = 4;
	for (int j = 0; j < node2.M; j++) {
		ASSERT_EQ(node(j), j);
	}
}

TEST(IdealElastic1DNode, Access)
{
	IdealElastic1DNode node;
	for (int i = 0; i < node.M; i++) {
		node(i) = i;
	}
	ASSERT_EQ(node.Vx, 0);
	ASSERT_EQ(node.Sxx, 1);

	IdealElastic1DNode node2;
	node2.Vx = 0; node2.Sxx = 1;
	for (int j = 0; j < node2.M; j++) {
		ASSERT_EQ(node(j), j);
	}
}