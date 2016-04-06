// #include <gtest/gtest.h>

// #include <lib/util/Area.hpp>

// #include <test/wrappers/Wrappers.hpp>
// #include <lib/rheology/models/Model.hpp>

// using namespace gcm;

// TEST(Solver, StageXForward)
// {
//	for (int accuracyOrder = 1; accuracyOrder < 5; accuracyOrder++) {
//		Task task;
//		Statement statement;
//		task.cubicGrid.dimensionality = 2;
//		task.cubicGrid.borderSize = accuracyOrder;
//		statement.globalSettings.CourantNumber = 1.0;
//		statement.materialConditions.defaultMaterial =
// std::make_shared<IsotropicMaterial>(4, 2, 0.5);
//		task.cubicGrid.sizes = {10, 10, 1};
//		task.cubicGrid.lengths = {2, 3, 1};
//		statement.globalSettings.numberOfSnaps = 1;

//		Statement::InitialCondition::Wave wave;
//		wave.waveType = Waves::T::P_FORWARD;
//		wave.direction = 0; // along x
//		wave.quantity = PhysicalQuantities::T::PRESSURE;
//		wave.quantityValue = 5;
//		Real3 min({0.3, -1, -1});
//		Real3 max({0.7, 4, 1});
//		wave.area = std::make_shared<AxisAlignedBoxArea>(min, max);
//		statement.initialCondition.waves.push_back(wave);

//		task.statements.push_back(statement);

//		CubicGrid::preprocessTask(task.cubicGrid);
//		DefaultSolverWrapper<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>>
// solver(task);
//		solver.beforeStatement(statement);
//		DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>::PdeVector zero({0, 0, 0,
// 0, 0});
//		auto pWave = solver.getMesh()->pde({2, 0, 0});

//		for (int i = 0; i < 7; i++) {
//			for (int y = 0; y < task.cubicGrid.sizes(1); y++) {
//				for (int x = 0; x < task.cubicGrid.sizes(0); x++) {
//					ASSERT_EQ(solver.getMesh()->pde({x, y, 0}),
//					           (x == 2 + i || x == 3 + i) ? pWave : zero)
//					<< "accuracyOrder = " << accuracyOrder << " i = " << i << "
// y = " << y << " x = " << x;
//				}
//			}
//			solver.stageForTest(0, solver.calculateTimeStep());
//		}
//	}
// }


// TEST(Solver, StageY)
// {
//	for (int accuracyOrder = 1; accuracyOrder < 5; accuracyOrder++) {
//		Task task;
//		Statement statement;
//		task.cubicGrid.dimensionality = 2;
//		task.cubicGrid.borderSize = accuracyOrder;
//		statement.globalSettings.CourantNumber = 1.0;
//		statement.materialConditions.defaultMaterial =
// std::make_shared<IsotropicMaterial>(4, 2, 0.5);
//		task.cubicGrid.sizes = {10, 10, 1};
//		task.cubicGrid.lengths = {3, 2, 1};
//		statement.globalSettings.numberOfSnaps = 1;

//		Statement::InitialCondition::Wave wave;
//		wave.waveType = Waves::T::P_FORWARD;
//		wave.direction = 1; // along y
//		wave.quantity = PhysicalQuantities::T::Vy;
//		wave.quantityValue = -2;
//		Real3 min({ -1, 0.3, -1});
//		Real3 max({ 4, 0.7, 1});
//		wave.area = std::make_shared<AxisAlignedBoxArea>(min, max);
//		statement.initialCondition.waves.push_back(wave);

//		task.statements.push_back(statement);

//		CubicGrid::preprocessTask(task.cubicGrid);
//		DefaultSolverWrapper<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>>
// solver(task);
//		solver.beforeStatement(statement);
//		auto pWave = solver.getMesh()->pde({0, 2, 0});
//		DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>::PdeVector zero({0, 0, 0,
// 0, 0});

//		for (int i = 0; i < 2; i++) {
//			for (int y = 0; y < task.cubicGrid.sizes(1); y++) {
//				for (int x = 0; x < task.cubicGrid.sizes(0); x++) {
//					ASSERT_EQ(solver.getMesh()->pde({x, y, 0}),
//					          (y == 2 + i || y == 3 + i) ? pWave : zero)
//					<< "accuracyOrder = " << accuracyOrder << " i = " << i << "
// y = " << y << " x = " << x;
//				}
//			}
//			solver.stageForTest(1, solver.calculateTimeStep());
//		}
//	}
// }


// TEST(Solver, StageYSxx)
// {
//	for (int accuracyOrder = 1; accuracyOrder < 5; accuracyOrder++) {
//		Task task;
//		Statement statement;
//		task.cubicGrid.dimensionality = 2;
//		task.cubicGrid.borderSize = accuracyOrder;
//		statement.globalSettings.CourantNumber = 0.7;
//		statement.materialConditions.defaultMaterial =
// std::make_shared<IsotropicMaterial>(4, 2, 0.5);
//		task.cubicGrid.sizes = {11, 11, 1};
//		task.cubicGrid.lengths = {1, 1, 1};
//		statement.globalSettings.numberOfSnaps = 1;

//		Statement::InitialCondition::Quantity quantity;
//		quantity.physicalQuantity = PhysicalQuantities::T::Sxx;
//		quantity.value = 10;
//		Real3 begin({0.5, 0.5, -1});
//		Real3 end({0.5, 0.5, 1});
//		quantity.area = std::make_shared<StraightBoundedCylinderArea>(0.01, begin, end);
//		statement.initialCondition.quantities.push_back(quantity);

//		task.statements.push_back(statement);

//		CubicGrid::preprocessTask(task.cubicGrid);
//		DefaultSolverWrapper<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>>
// solver(task);
//		solver.beforeStatement(statement);
//		auto sxxOnly = solver.getMesh()->pde({5, 5, 0});
//		DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>::PdeVector zero({0, 0, 0,
// 0, 0});

//		for (int i = 0; i < 7; i++) {
//			for (int y = 0; y < task.cubicGrid.sizes(1); y++) {
//				for (int x = 0; x < task.cubicGrid.sizes(0); x++) {
//					ASSERT_EQ(solver.getMesh()->pde({x, y, 0}), (x == 5 && y ==
// 5 ) ? sxxOnly : zero)
//					<< "accuracyOrder = " << accuracyOrder << " i = " << i << "
// y = " << y << " x = " << x;
//				}
//			}
//			solver.stageForTest(1, solver.calculateTimeStep());
//		}
//	}
// }
