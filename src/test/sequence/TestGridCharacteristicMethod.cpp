#include <gtest/gtest.h>

#include <lib/util/Area.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/mesh/DefaultMesh.hpp>
#include <lib/numeric/gcm/GridCharacteristicMethodCubicGrid.hpp>
#include <lib/rheology/models/models.hpp>

#include <lib/util/snapshot/VtkSnapshotter.hpp>

using namespace gcm;


TEST(GridCharacteristicMethodCubicGrid, interpolateValuesAround) {
	Task task;
	
	task.materialConditions.byAreas.defaultMaterial =
			std::make_shared<IsotropicMaterial>(2, 2, 1);
	
	task.cubicGrid.sizes = {3, 3};
	task.cubicGrid.borderSize = 1;
	task.cubicGrid.lengths = {2, 2}; // h_x = h_y = 1.0
	
	Task::InitialCondition::Quantity quantity;
	quantity.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	quantity.value = -1.0;
	quantity.area = std::make_shared<SphereArea>(0.1, Real3({1, 1, 0}));
	task.initialCondition.quantities.push_back(quantity);
	
	
	for (int stage = 0; stage <= 1; stage++) {
		
		typedef CubicGrid<2> Grid;
		typedef typename Grid::GlobalScene GS;
		std::shared_ptr<GS> gs(new GS(task));
		
		DefaultMesh<ElasticModel<2>, Grid, IsotropicMaterial> mesh(task, gs.get(), 0);
		for (int x = 0; x < task.cubicGrid.sizes.at(0); x++) {
			for (int y = 0; y < task.cubicGrid.sizes.at(1); y++) {
				// check that values is set properly
				ASSERT_EQ(mesh.pdeVars({x, y}).velocity(0), 0.0);
				ASSERT_EQ(mesh.pdeVars({x, y}).velocity(1), 0.0);
				ASSERT_EQ(mesh.pdeVars({x, y}).sigma(0, 0),
				          (x == 1 && y == 1) ? 1.0 : 0.0);
				ASSERT_EQ(mesh.pdeVars({x, y}).sigma(0, 1), 0.0);
				ASSERT_EQ(mesh.pdeVars({x, y}).sigma(1, 1),
				          (x == 1 && y == 1) ? 1.0 : 0.0);
			}
		}
		
		auto m = GridCharacteristicMethod<
				ElasticModel<2>, Grid, IsotropicMaterial>().
						interpolateValuesAround(
								mesh, stage, {1, 1}, {-1, 1, -0.5, 0.5, 0});
		
		for (int i = 0; i < 5; i++) {
			ASSERT_EQ(m(i, 0), 0.0) << "i = " << i; // Courant = 1
			ASSERT_EQ(m(i, 1), 0.0) << "i = " << i; // Courant = 1
			ASSERT_EQ(m(0, i), 0.0) << "i = " << i; // Vx
			ASSERT_EQ(m(1, i), 0.0) << "i = " << i; // Vy
			ASSERT_EQ(m(3, i), 0.0) << "i = " << i; // Sxy
		}
		
		ASSERT_EQ(m(2, 2), 0.5); // Courant = 0.5
		ASSERT_EQ(m(2, 3), 0.5); // Courant = 0.5
		ASSERT_EQ(m(4, 2), 0.5); // Courant = 0.5
		ASSERT_EQ(m(4, 3), 0.5); // Courant = 0.5
		
		ASSERT_EQ(m(2, 4), 1.0); // Courant = 0
		ASSERT_EQ(m(4, 4), 1.0); // Courant = 0
	}
}


TEST(GridCharacteristicMethodCubicGrid, StageYForward) {
	for (int accuracyOrder = 1; accuracyOrder < 5; accuracyOrder++) {
		Task task;
		
		task.bodies = {{0, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}}};
		task.cubicGrid.borderSize = accuracyOrder;
		task.globalSettings.CourantNumber = 1;
		task.materialConditions.byAreas.defaultMaterial = 
				std::make_shared<IsotropicMaterial>(4, 2, 0.5);
		task.cubicGrid.sizes = {10, 10};
		task.cubicGrid.lengths = {3, 2};
	
		Task::InitialCondition::Wave wave;
		wave.waveType = Waves::T::P_FORWARD;
		wave.direction = 1; // along y
		wave.quantity = PhysicalQuantities::T::PRESSURE;
		wave.quantityValue = 5;
		Real3 min({ -1, 0.3, -1});
		Real3 max({  4, 0.7, 1});
		wave.area = std::make_shared<AxisAlignedBoxArea>(min, max);
		task.initialCondition.waves.push_back(wave);
		
		typedef CubicGrid<2> Grid;
		typedef typename Grid::GlobalScene GS;
		std::shared_ptr<GS> gs(new GS(task));
		
		DefaultSolver<DefaultMesh<
				ElasticModel<2>, CubicGrid<2>, IsotropicMaterial>>
						solver(task, gs.get(), 0);
		
		const auto zero = ElasticModel<2>::PdeVector::Zeros();
		auto pWave = solver.getMesh()->pde({0, 2});
		
		for (int i = 0; i < 7; i++) {
			for (int y = 0; y < task.cubicGrid.sizes.at(1); y++) {
				for (int x = 0; x < task.cubicGrid.sizes.at(0); x++) {
					auto value = solver.getMesh()->pde({x, y});
					auto analytic = (y == 2 + i || y == 3 + i) ? pWave : zero;
					ASSERT_TRUE(linal::approximatelyEqual(value, analytic))
						<< "accuracyOrder = " << accuracyOrder << " i = " << i 
						<< " y = " << y << " x = " << x
						<< "value" << value;
				}
			}
			solver.privateStage(1, solver.calculateTimeStep());
		}
	}
}




