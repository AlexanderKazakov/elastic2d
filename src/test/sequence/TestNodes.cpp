#include <gtest/gtest.h>

#include "lib/nodes/IdealElastic3DNode.hpp"
#include "lib/nodes/IdealElastic2DNode.hpp"
#include "lib/nodes/IdealElastic1DNode.hpp"

using namespace gcm;

TEST(IdealElastic3DNode, Access)
{
	IdealElastic3DNode node;
	for (int i = 0; i < node.u.M; i++) {
		node.u(i) = i;
	}
	ASSERT_EQ(node.u.Vx, 0);
	ASSERT_EQ(node.u.Vy, 1);
	ASSERT_EQ(node.u.Vz, 2);
	ASSERT_EQ(node.u.Sxx, 3);
	ASSERT_EQ(node.u.Sxy, 4);
	ASSERT_EQ(node.u.Sxz, 5);
	ASSERT_EQ(node.u.Syy, 6);
	ASSERT_EQ(node.u.Syz, 7);
	ASSERT_EQ(node.u.Szz, 8);

	IdealElastic3DNode node2;
	node2.u.Vx = 0; node2.u.Vy = 1; node2.u.Vz = 2;
	node2.u.Sxx = 3; node2.u.Sxy = 4; node2.u.Sxz = 5; node2.u.Syy = 6; node2.u.Syz = 7; node2.u.Szz = 8;
	for (int j = 0; j < node2.u.M; j++) {
		ASSERT_EQ(node.u(j), j);
	}
}

TEST(IdealElastic2DNode, Access)
{
	IdealElastic2DNode node;
	for (int i = 0; i < node.u.M; i++) {
		node.u(i) = i;
	}
	ASSERT_EQ(node.u.Vx, 0);
	ASSERT_EQ(node.u.Vy, 1);
	ASSERT_EQ(node.u.Sxx, 2);
	ASSERT_EQ(node.u.Sxy, 3);
	ASSERT_EQ(node.u.Syy, 4);

	IdealElastic2DNode node2;
	node2.u.Vx = 0; node2.u.Vy = 1;
	node2.u.Sxx = 2; node2.u.Sxy = 3; node2.u.Syy = 4;
	for (int j = 0; j < node2.u.M; j++) {
		ASSERT_EQ(node.u(j), j);
	}
}

TEST(IdealElastic1DNode, Access)
{
	IdealElastic1DNode node;
	for (int i = 0; i < node.u.M; i++) {
		node.u(i) = i;
	}
	ASSERT_EQ(node.u.Vx, 0);
	ASSERT_EQ(node.u.Sxx, 1);

	IdealElastic1DNode node2;
	node2.u.Vx = 0; node2.u.Sxx = 1;
	for (int j = 0; j < node2.u.M; j++) {
		ASSERT_EQ(node.u(j), j);
	}
}