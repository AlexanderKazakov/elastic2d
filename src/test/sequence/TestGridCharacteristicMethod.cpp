#include <gtest/gtest.h>

#include <libgcm/util/math/Area.hpp>
#include <libgcm/grid/cubic/CubicGrid.hpp>
#include <libgcm/engine/mesh/DefaultMesh.hpp>
#include <libgcm/engine/cubic/GridCharacteristicMethod.hpp>
#include <libgcm/rheology/models/models.hpp>

#include <libgcm/util/snapshot/VtkSnapshotter.hpp>

using namespace gcm;


TEST(GridCharacteristicMethodCubicGrid, interpolateValuesAround) {
	Task task;
	task.materialConditions.byAreas.defaultMaterial =
			std::make_shared<IsotropicMaterial>(2, 2, 1);
	
	Task::InitialCondition::Quantity quantity;
	quantity.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	quantity.value = -1.0;
	quantity.area = std::make_shared<SphereArea>(0.1, Real3({1, 1, 0}));
	task.initialCondition.quantities.push_back(quantity);
	
	typedef CubicGrid<2> Grid;
	typedef typename Grid::ConstructionPack ConstructionPack;
	ConstructionPack cp;
	cp.borderSize = 1;
	cp.sizes = {3, 3};
	cp.h = {1, 1};
	
	for (int stage = 0; stage <= 1; stage++) {
		
		DefaultMesh<ElasticModel<2>, Grid, IsotropicMaterial> mesh(task, 0, cp, 1);
		for (int x = 0; x < cp.sizes(0); x++) {
			for (int y = 0; y < cp.sizes(1); y++) {
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
		
		auto m = cubic::GridCharacteristicMethod<
				DefaultMesh<ElasticModel<2>, Grid, IsotropicMaterial>>(task).
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


