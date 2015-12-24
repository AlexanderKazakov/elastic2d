#include <gtest/gtest.h>

#include "lib/MPISolver.hpp"

using namespace gcm;

TEST(Solver, StageXForward)
{
	for (int accuracyOrder = 1; accuracyOrder < 20; accuracyOrder++) {
		Task task;
		task.accuracyOrder = accuracyOrder;
		task.CourantNumber = 1.0;
		task.lambda0 = 2.0;
		task.mu0 = 0.5;
		task.rho0 = 4.0;
		task.X = 10;
		task.Y = 10;
		task.xLength = 2.0;
		task.yLength = 3.0;
		task.numberOfSnaps = 1;
		task.T = 100.0;
		task.initialConditions = InitialConditions::PWaveX;

		Mesh *mesh = new Mesh();
		Mesh *newMesh = new Mesh();
		mesh->initialize(task);
		newMesh->initialize(task);
		Vector pWave = mesh->getNodeForTest(0, 2).u;
		Vector zero; zero.createVector({0, 0, 0, 0, 0});
		MPISolver solver(mesh, newMesh);
		for (int i = 0; i < 7; i++) {
			for (int y = 0; y < task.Y; y++) {
				for (int x = 0; x < task.X; x++) {
					ASSERT_EQ(mesh->getNodeForTest(y, x).u,
					          (x == 2 + i || x == 3 + i) ? pWave : zero)
					<< "accuracyOrder = " << accuracyOrder << " i = " << i << " y = " << y << " x = " << x;
				}
			}
			solver.stage(0, mesh->getTauForTest());
			std::swap(mesh, newMesh); // because solver swap them internally
		}
	}
}


TEST(Solver, StageXBackward)
{
/*
	for (int accuracyOrder = 2; accuracyOrder < 20; accuracyOrder++) {
		Task task;
		task.accuracyOrder = accuracyOrder;
		task.lambda0 = 2.0;
		task.mu0 = task.lambda0 / 2; // s-wave two times slower than p-wave
		task.CourantNumber = 2.0; // so Courant for s-wave is 1.0
		task.rho0 = 4.0;
		task.X = 10;
		task.Y = 10;
		task.xLength = 2.0;
		task.yLength = 3.0;
		task.numberOfSnaps = 1;
		task.T = 100.0;
		task.initialConditions = InitialConditions::SWaveXBackward;

		Mesh *mesh = new Mesh();
		Mesh *newMesh = new Mesh();
		mesh->initialize(task);
		newMesh->initialize(task);
		Vector pWave = mesh->getNodeForTest(0, task.X - 3).u;
		Vector zero; zero.createVector({0, 0, 0, 0, 0});
		MPISolver solver(mesh, newMesh);
		for (int i = 0; i < 7; i++) {
			for (int y = 0; y < task.Y; y++) {
				for (int x = 0; x < task.X; x++) {
					ASSERT_EQ(mesh->getNodeForTest(y, x).u,
					          (x == task.X - 4 - i || x == task.X - 3 - i) ? pWave : zero)
					<< "accuracyOrder = " << accuracyOrder << " i = " << i << " y = " << y << " x = " << x;
				}
			}
			solver.stage(0, mesh->getTauForTest());
			std::swap(mesh, newMesh); // because solver swap them internally
		}
	}
*/
}


TEST(Solver, StageY)
{
	for (int accuracyOrder = 1; accuracyOrder < 20; accuracyOrder++) {
		Task task;
		task.accuracyOrder = accuracyOrder;
		task.CourantNumber = 1.0;
		task.lambda0 = 2.0;
		task.mu0 = 0.5;
		task.rho0 = 4.0;
		task.X = 10;
		task.Y = 10;
		task.xLength = 3.0;
		task.yLength = 3.0;
		task.numberOfSnaps = 1;
		task.T = 100.0;
		task.initialConditions = InitialConditions::PWaveY;

		Mesh *mesh = new Mesh();
		Mesh *newMesh = new Mesh();
		mesh->initialize(task);
		newMesh->initialize(task);
		Vector pWave = mesh->getNodeForTest(2, 0).u;
		Vector zero; zero.createVector({0, 0, 0, 0, 0});
		MPISolver solver(mesh, newMesh);
		for (int i = 0; i < 2; i++) {
			for (int y = 0; y < task.Y; y++) {
				for (int x = 0; x < task.X; x++) {
					ASSERT_EQ(mesh->getNodeForTest(y, x).u,
					          (y == 2 + i || y == 3 + i || y == 4 + i || y == 5 + i || y == 6 + i) ? pWave : zero)
					<< "accuracyOrder = " << accuracyOrder << " i = " << i << " y = " << y << " x = " << x;
				}
			}
			solver.stage(1, mesh->getTauForTest());
			std::swap(mesh, newMesh); // because solver swap them internally
		}
	}
}


TEST(Solver, StageYSxx)
{
	for (int accuracyOrder = 1; accuracyOrder < 20; accuracyOrder++) {
		Task task;
		task.accuracyOrder = accuracyOrder;
		task.CourantNumber = 0.7;
		task.lambda0 = 2.0;
		task.mu0 = 0.5;
		task.rho0 = 4.0;
		task.X = 20;
		task.Y = 10;
		task.xLength = 7.0;
		task.yLength = 3.0;
		task.numberOfSnaps = 1;
		task.T = 100.0;
		task.initialConditions = InitialConditions::SxxOnly;

		Mesh *mesh = new Mesh();
		Mesh *newMesh = new Mesh();
		mesh->initialize(task);
		newMesh->initialize(task);
		Vector sxxOnly = mesh->getNodeForTest(task.Y / 2, task.X / 2).u;
		Vector zero; zero.createVector({0, 0, 0, 0, 0});
		MPISolver solver(mesh, newMesh);
		for (int i = 0; i < 7; i++) {
			for (int y = 0; y < task.Y; y++) {
				for (int x = 0; x < task.X; x++) {
					ASSERT_EQ(mesh->getNodeForTest(y, x).u,
					          (x == task.X / 2 && y == task.Y / 2 ) ? sxxOnly : zero)
					<< "accuracyOrder = " << accuracyOrder << " i = " << i << " y = " << y << " x = " << x;
				}
			}
			solver.stage(1, mesh->getTauForTest());
			std::swap(mesh, newMesh); // because solver swap them internally
		}
	}
}


TEST(Solver, calculate)
{
	Task task;
	task.accuracyOrder = 5;
	task.CourantNumber = 4.5;
	task.lambda0 = 2.0;
	task.mu0 = 0.5;
	task.rho0 = 4.0;
	task.X = 20;
	task.Y = 40;
	task.xLength = 7.0;
	task.yLength = 3.0;
	task.numberOfSnaps = 9;
	task.T = 100.0;
	task.initialConditions = InitialConditions::SWaveY;

	Mesh *mesh = new Mesh();
	Mesh *newMesh = new Mesh();
	mesh->initialize(task);
	newMesh->initialize(task);
	Vector sWave = mesh->getNodeForTest(3, task.X / 2).u;
	MPISolver solver(mesh, newMesh);
	solver.calculate();
	ASSERT_EQ(sWave, mesh->getNodeForTest(22, task.X / 2).u);
}


TEST(Solver, TwoLayersDifferentRho)
{
	real rho2rho0Initial = 0.25;
	int numberOfSnapsInitial = 30;
	for (int i = 0; i < 5; i++) {

		Task task;
		task.accuracyOrder = 3;
		task.CourantNumber = 1.5;
		task.lambda0 = 2.0;
		task.mu0 = 0.8;
		task.rho0 = 1.0;
		task.X = 50;
		task.Y = 100;
		task.xLength = 2.0;
		task.yLength = 1.0;
		task.numberOfSnaps = numberOfSnapsInitial + 2 * i; // in order to catch the impulses
		task.initialConditions = InitialConditions::PWaveY;

		Mesh *mesh = new Mesh();
		Mesh *newMesh = new Mesh();
		mesh->initialize(task);
		newMesh->initialize(task);

		real rho2rho0 = rho2rho0Initial * pow(2, i);
		real lambda2lambda0 = 1;
		real mu2mu0 = 1;
		mesh->changeRheology(rho2rho0, lambda2lambda0, mu2mu0);
		newMesh->changeRheology(rho2rho0, lambda2lambda0, mu2mu0);


		int leftNodeIndex = task.Y * 0.25;
		Node init = mesh->getNodeForTest(leftNodeIndex, task.X / 2);

		MPISolver solver(mesh, newMesh);
		solver.calculate();

		int rightNodeIndex = task.Y * 0.7;
		Node reflect = mesh->getNodeForTest(leftNodeIndex, task.X / 2);
		Node transfer = mesh->getNodeForTest(rightNodeIndex, task.X / 2);

		auto defaultMatrix = mesh->getDefaultMatrixForTest();
		real rho0 = defaultMatrix->rho;
		real lambda0 = defaultMatrix->lambda;
		real mu0 = defaultMatrix->mu;
		real E0 = mu0 * (3 * lambda0 + 2 * mu0) / (lambda0 + mu0); // Young's modulus
		real Z0 = sqrt(E0 * rho0); // acoustic impedance


		real rho = rho2rho0 * rho0;
		real lambda = lambda2lambda0 * lambda0;
		real mu = mu2mu0 * mu0;
		real E = mu * (3 * lambda + 2 * mu) / (lambda + mu); // Young's modulus
		real Z = sqrt(E * rho); // acoustic impedance

		ASSERT_NEAR(reflect.get(NodeMap::Syy) / init.get(NodeMap::Syy),
		            (Z - Z0) / (Z + Z0),
		            1e-2);
		ASSERT_NEAR(reflect.get(NodeMap::Vy) / init.get(NodeMap::Vy),
		            (Z0 - Z) / (Z + Z0),
		            1e-2);

		ASSERT_NEAR(transfer.get(NodeMap::Syy) / init.get(NodeMap::Syy),
		            2 * Z / (Z + Z0),
		            1e-2);
		ASSERT_NEAR(transfer.get(NodeMap::Vy) / init.get(NodeMap::Vy),
		            2 * Z0 / (Z + Z0),
		            1e-2);


	}
}


TEST(Solver, TwoLayersDifferentE)
{
	real E2E0Initial = 0.25;
	int numberOfSnapsInitial = 40;
	for (int i = 0; i < 5; i++) {

		Task task;
		task.accuracyOrder = 3;
		task.CourantNumber = 1.5;
		task.lambda0 = 2.0;
		task.mu0 = 0.8;
		task.rho0 = 1.0;
		task.X = 50;
		task.Y = 100;
		task.xLength = 2.0;
		task.yLength = 1.0;
		task.numberOfSnaps = numberOfSnapsInitial - 2 * i; // in order to catch the impulses
		task.initialConditions = InitialConditions::PWaveY;

		Mesh *mesh = new Mesh();
		Mesh *newMesh = new Mesh();
		mesh->initialize(task);
		newMesh->initialize(task);

		real rho2rho0 = 1;
		real lambda2lambda0 = E2E0Initial * pow(2, i);
		real mu2mu0 = E2E0Initial * pow(2, i);
		mesh->changeRheology(rho2rho0, lambda2lambda0, mu2mu0);
		newMesh->changeRheology(rho2rho0, lambda2lambda0, mu2mu0);


		int leftNodeIndex = task.Y * 0.25;
		Node init = mesh->getNodeForTest(leftNodeIndex, task.X / 2);

		MPISolver solver(mesh, newMesh);
		solver.calculate();

		int rightNodeIndex = task.Y * 0.7;
		Node reflect = mesh->getNodeForTest(leftNodeIndex, task.X / 2);
		Node transfer = mesh->getNodeForTest(rightNodeIndex, task.X / 2);

		auto defaultMatrix = mesh->getDefaultMatrixForTest();
		real rho0 = defaultMatrix->rho;
		real lambda0 = defaultMatrix->lambda;
		real mu0 = defaultMatrix->mu;
		real E0 = mu0 * (3 * lambda0 + 2 * mu0) / (lambda0 + mu0); // Young's modulus
		real Z0 = sqrt(E0 * rho0); // acoustic impedance


		real rho = rho2rho0 * rho0;
		real lambda = lambda2lambda0 * lambda0;
		real mu = mu2mu0 * mu0;
		real E = mu * (3 * lambda + 2 * mu) / (lambda + mu); // Young's modulus
		real Z = sqrt(E * rho); // acoustic impedance

		ASSERT_NEAR(reflect.get(NodeMap::Syy) / init.get(NodeMap::Syy),
		            (Z - Z0) / (Z + Z0),
		            1e-2);
		ASSERT_NEAR(reflect.get(NodeMap::Vy) / init.get(NodeMap::Vy),
		            (Z0 - Z) / (Z + Z0),
		            1e-2);

		ASSERT_NEAR(transfer.get(NodeMap::Syy) / init.get(NodeMap::Syy),
		            2 * Z / (Z + Z0),
		            1e-2);
		ASSERT_NEAR(transfer.get(NodeMap::Vy) / init.get(NodeMap::Vy),
		            2 * Z0 / (Z + Z0),
		            1e-2);


	}
}