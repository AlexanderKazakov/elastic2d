#include <gtest/gtest.h>

#include <lib/util/areas/AxisAlignedBoxArea.hpp>
#include <lib/util/areas/StraightBoundedCylinderArea.hpp>

#include <test/wrappers/Wrappers.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

TEST(Solver, StageXForward)
{
	for (int accuracyOrder = 1; accuracyOrder < 5; accuracyOrder++) {
		Task task;
		task.accuracyOrder = accuracyOrder;
		task.CourantNumber = 1.0;
		task.isotropicMaterial = IsotropicMaterial(4.0, 2.0, 0.5);
		task.sizes(0) = 10;
		task.sizes(1) = 10;
		task.lengthes = {2, 3, 1};
		task.numberOfSnaps = 1;
		task.T = 100.0;

		Task::InitialCondition::Wave wave;
		wave.waveType = Waves::T::P_FORWARD;
		wave.direction = 0; // along x
		wave.quantity = PhysicalQuantities::T::PRESSURE;
		wave.quantityValue = 5;
		linal::Vector3 min({0.3, -1, -1});
		linal::Vector3 max({0.7, 4, 1});
		wave.area = std::make_shared<AxisAlignedBoxArea>(min, max);
		task.initialCondition.waves.push_back(wave);

		DefaultSolverWrapper<StructuredGrid<Elastic2DModel>> solver;
		solver.initialize(task);
		auto pWave = solver.getMesh()->pde(2, 0, 0);
		StructuredGrid<Elastic2DModel>::PdeVector zero({0, 0, 0, 0, 0});

		for (int i = 0; i < 7; i++) {
			for (int y = 0; y < task.sizes(1); y++) {
				for (int x = 0; x < task.sizes(0); x++) {
					ASSERT_EQ(solver.getMesh()->pde(x, y, 0),
					           (x == 2 + i || x == 3 + i) ? pWave : zero)
					<< "accuracyOrder = " << accuracyOrder << " i = " << i << " y = " << y << " x = " << x;
				}
			}
			solver.stageForTest(0, solver.getTauForTest());
		}
	}
}


TEST(Solver, StageY)
{
	for (int accuracyOrder = 1; accuracyOrder < 5; accuracyOrder++) {
		Task task;
		task.accuracyOrder = accuracyOrder;
		task.CourantNumber = 1.0;
		task.isotropicMaterial = IsotropicMaterial(4.0, 2.0, 0.5);
		task.sizes(0) = 10;
		task.sizes(1) = 10;
		task.lengthes = {3, 2, 1};
		task.numberOfSnaps = 1;
		task.T = 100.0;

		Task::InitialCondition::Wave wave;
		wave.waveType = Waves::T::P_FORWARD;
		wave.direction = 1; // along y
		wave.quantity = PhysicalQuantities::T::Vy;
		wave.quantityValue = -2;
		linal::Vector3 min({ -1, 0.3, -1});
		linal::Vector3 max({ 4, 0.7, 1});
		wave.area = std::make_shared<AxisAlignedBoxArea>(min, max);
		task.initialCondition.waves.push_back(wave);

		DefaultSolverWrapper<StructuredGrid<Elastic2DModel>> solver;
		solver.initialize(task);
		auto pWave = solver.getMesh()->pde(0, 2, 0);
		StructuredGrid<Elastic2DModel>::PdeVector zero({0, 0, 0, 0, 0});

		for (int i = 0; i < 2; i++) {
			for (int y = 0; y < task.sizes(1); y++) {
				for (int x = 0; x < task.sizes(0); x++) {
					ASSERT_EQ(solver.getMesh()->pde(x, y, 0),
					          (y == 2 + i || y == 3 + i) ? pWave : zero)
					<< "accuracyOrder = " << accuracyOrder << " i = " << i << " y = " << y << " x = " << x;
				}
			}
			solver.stageForTest(1, solver.getTauForTest());
		}
	}
}


TEST(Solver, StageYSxx)
{
	for (int accuracyOrder = 1; accuracyOrder < 5; accuracyOrder++) {
		Task task;
		task.accuracyOrder = accuracyOrder;
		task.CourantNumber = 0.7;
		task.isotropicMaterial = IsotropicMaterial(4.0, 2.0, 0.5);
		task.sizes(0) = 11;
		task.sizes(1) = 11;
		task.lengthes = {1, 1, 1};
		task.numberOfSnaps = 1;

		Task::InitialCondition::Quantity quantity;
		quantity.physicalQuantity = PhysicalQuantities::T::Sxx;
		quantity.value = 10;
		linal::Vector3 begin({0.5, 0.5, -1});
		linal::Vector3 end({0.5, 0.5, 1});
		quantity.area = std::make_shared<StraightBoundedCylinderArea>(0.01, begin, end);
		task.initialCondition.quantities.push_back(quantity);

		DefaultSolverWrapper<StructuredGrid<Elastic2DModel>> solver;
		solver.initialize(task);
		auto sxxOnly = solver.getMesh()->pde(5, 5, 0);
		StructuredGrid<Elastic2DModel>::PdeVector zero({0, 0, 0, 0, 0});

		for (int i = 0; i < 7; i++) {
			for (int y = 0; y < task.sizes(1); y++) {
				for (int x = 0; x < task.sizes(0); x++) {
					ASSERT_EQ(solver.getMesh()->pde(x, y, 0), (x == 5 && y == 5 ) ? sxxOnly : zero)
					<< "accuracyOrder = " << accuracyOrder << " i = " << i << " y = " << y << " x = " << x;
				}
			}
			solver.stageForTest(1, solver.getTauForTest());
		}
	}
}