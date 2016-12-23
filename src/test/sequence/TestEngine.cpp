#include <gtest/gtest.h>

#include <libgcm/engine/cubic/Engine.hpp>
#include <libgcm/util/math/Area.hpp>
#include <libgcm/rheology/models/models.hpp>

#include <libgcm/engine/mesh/DefaultMesh.hpp>
#include <libgcm/grid/cubic/CubicGrid.hpp>


using namespace gcm;
using namespace gcm::cubic;


struct Wrapper {
	typedef DefaultMesh<ElasticModel<2>, CubicGrid<2>, IsotropicMaterial> Mesh;
	static std::shared_ptr<const Mesh> getMesh(
			const Engine<2>& engine, const size_t id = 0) {
		auto grid = engine.getMesh(id);
		auto mesh = std::dynamic_pointer_cast<const Mesh>(grid);
		assert_true(mesh);
		return mesh;
	}
};


TEST(Engine, AdhesionContact) {
	int Y = 41, X = 21;
	
	Task task;
	task.globalSettings.dimensionality = 2;
	task.globalSettings.gridId = Grids::T::CUBIC;
	task.globalSettings.snapshottersId = { Snapshotters::T::VTK };
	task.globalSettings.numberOfSnaps = 70;
	task.globalSettings.CourantNumber = 0.9;
	
	task.bodies = {
			{0, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}},
			{1, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}}
	};
	
	task.materialConditions.byAreas.defaultMaterial = 
			std::make_shared<IsotropicMaterial>(4, 2, 0.5);
	task.cubicGrid.borderSize = 2;
	task.cubicGrid.h = {1, 0.25};
	task.cubicGrid.cubics = {
			{0, {{X, Y}, { 0, 0}}},
			{1, {{X, Y}, { 0, Y}}}
	};
	
	Task::InitialCondition::Wave wave;
	wave.waveType = Waves::T::P_FORWARD;
	wave.direction = 1; // along y
	wave.quantity = PhysicalQuantities::T::PRESSURE;
	wave.quantityValue = 1;
	Real3 min({-1000, 2.5, -1000});
	Real3 max({ 1000, 7.5,  1000});
	wave.area = std::make_shared<AxisAlignedBoxArea>(min, max);
	task.initialCondition.waves.push_back(wave);
	
	Engine<2> two(task);
	two.run();
	auto first  = Wrapper::getMesh(two, 0);
	auto second = Wrapper::getMesh(two, 1);
	
	
	task.bodies = {
			{0, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}},
	};
	task.cubicGrid.cubics = {
			{0, {{X, 2 * Y}, { 0,  0}}},
	};
	Engine<2> one(task);
	one.run();
	auto all = Wrapper::getMesh(one, 0);
	
	
	for (int x = 0; x < X; x++) {
		for (int y = 0; y < Y; y++) {
			ASSERT_EQ(all->coords({x, y}), first->coords({x, y}));
			ASSERT_EQ(all->coords({x, Y + y}), second->coords({x, y}));
			ASSERT_EQ(all->pde({x, y}), first->pde({x, y}));
			ASSERT_EQ(all->pde({x, Y + y}), second->pde({x, y}));
		}
	}
}



TEST(Engine, runStatement) {
	Task task;
	task.globalSettings.dimensionality = 2;
	task.globalSettings.gridId = Grids::T::CUBIC;
	
	task.bodies = {
		{0, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}}
	};
	
	
	task.globalSettings.CourantNumber = 4.5;
	task.materialConditions.byAreas.defaultMaterial = 
			std::make_shared<IsotropicMaterial>(4, 2, 0.5);
	task.cubicGrid.borderSize = 5;
	task.cubicGrid.h = {7.0/19, 3.0/39};
	task.cubicGrid.cubics = {{0, {{20, 40}, {0, 0}}}};
	
	task.globalSettings.numberOfSnaps = 9;
	task.globalSettings.requiredTime = 100.0;
	
	Task::InitialCondition::Wave wave;
	wave.waveType = Waves::T::S1_FORWARD;
	wave.direction = 1; // along y
	wave.quantity = PhysicalQuantities::T::Vx;
	wave.quantityValue = 1;
	Real3 min({ -1, 0.1125, -1});
	Real3 max({ 8, 0.6375, 1});
	wave.area = std::make_shared<AxisAlignedBoxArea>(min, max);
	task.initialCondition.waves.push_back(wave);
	
//	try {
	Engine<2> engine(task);
	
	// s-wave
	auto expected = Wrapper::getMesh(engine)->pde({10, 3});
	ASSERT_NE(linal::zeros(expected), expected);
	engine.run();
	auto actual = Wrapper::getMesh(engine)->pde({10, 22});
	
	ASSERT_TRUE(linal::approximatelyEqual(expected, actual)) << expected << actual;
	
//	} catch (Exception e) {
//		std::cout << e.what();
//	}
}


TEST(Engine, TwoLayersDifferentRho) {
	real rho2rho0Initial = 0.25;
	for (int i = 0; i < 5; i++) {
		
		Task task;
		task.globalSettings.dimensionality = 2;
		task.globalSettings.gridId = Grids::T::CUBIC;
		
		task.bodies = {
			{0, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}}
		};
		
		task.globalSettings.CourantNumber = 1.5;
		task.globalSettings.numberOfSnaps = 0;
		task.globalSettings.requiredTime = 0.24;
		
		real rho0 = 1, lambda0 = 2, mu0 = 0.8;
		real rho2rho0 = rho2rho0Initial * pow(2, i), lambda2lambda0 = 1, mu2mu0 = 1;
		real rho = rho2rho0 * rho0, lambda = lambda2lambda0 * lambda0, mu = mu2mu0 * mu0;
		
		task.materialConditions.byAreas.defaultMaterial = 
				std::make_shared<IsotropicMaterial>(rho0, lambda0, mu0);
		Task::MaterialCondition::ByAreas::Inhomogenity newMaterial;
		newMaterial.area = std::make_shared<AxisAlignedBoxArea>(
				Real3({-10, 0.5 - 1e-5,-10}), Real3({10, 10, 10}));
		newMaterial.material = std::make_shared<IsotropicMaterial>(rho, lambda, mu);
		task.materialConditions.byAreas.materials.push_back(newMaterial);
		
		task.cubicGrid.borderSize = 3;
		task.cubicGrid.h = {2.0/49, 1.0/99};
		task.cubicGrid.cubics = {{0, {{50, 100}, {0, 0}}}};
		
		
		Task::InitialCondition::Wave wave;
		wave.waveType = Waves::T::P_FORWARD;
		wave.direction = 1; // along y
		wave.quantity = PhysicalQuantities::T::Vy;
		wave.quantityValue = -2;
		Real3 min({ -1, 0.015, -1});
		Real3 max({ 4, 0.455, 1});
		wave.area = std::make_shared<AxisAlignedBoxArea>(min, max);
		task.initialCondition.waves.push_back(wave);
		
		
		Engine<2> engine(task);
		auto init = Wrapper::getMesh(engine)->pdeVars({25, 25});
		ASSERT_NE(linal::zeros(init), init);
		
		engine.run();
		
//		int rightNodeIndex = (int) (task.cubicGrid.sizes(1) * 0.7);
		auto reflect = Wrapper::getMesh(engine)->pdeVars({25, 25});
				
//		auto transfer = Wrapper::getMesh(engine)->pdeVars(
//				{task.cubicGrid.sizes(0) / 2, rightNodeIndex});
		
		real E0 = mu0 * (3 * lambda0 + 2 * mu0) / (lambda0 + mu0); // Young's modulus
		real Z0 = sqrt(E0 * rho0);                                 // acoustic impedance
		real E = mu * (3 * lambda + 2 * mu) / (lambda + mu);       // Young's modulus
		real Z = sqrt(E * rho);                                    // acoustic impedance
		
		ASSERT_NEAR(reflect.sigma(1, 1) / init.sigma(1, 1),
		            (Z - Z0) / (Z + Z0),
		            1e-2);
		ASSERT_NEAR(reflect.velocity(1) / init.velocity(1),
		            (Z0 - Z) / (Z + Z0),
		            1e-2);
//		TODO - complete stable beautiful tests
//		ASSERT_NEAR(transfer.sigma(1, 1) / init.sigma(1, 1),
//		            2 * Z / (Z + Z0),
//		            1e-2);
//		ASSERT_NEAR(transfer.velocity(1) / init.velocity(1),
//		            2 * Z0 / (Z + Z0),
//		            1e-2);
		
		
	}
}


TEST(Engine, TwoLayersDifferentE) {
	real E2E0Initial = 0.25;
	for (int i = 0; i < 5; i++) {
		
		Task task;
		task.globalSettings.dimensionality = 2;
		task.globalSettings.gridId = Grids::T::CUBIC;
		task.globalSettings.CourantNumber = 1.5;
		
		task.bodies = {
			{0, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}}
		};
		
		
		real rho0 = 1, lambda0 = 2, mu0 = 0.8;
		real rho2rho0 = 1, lambda2lambda0 = E2E0Initial * pow(2, i), mu2mu0 = E2E0Initial * pow(2, i);
		real rho = rho2rho0 * rho0, lambda = lambda2lambda0 * lambda0, mu = mu2mu0 * mu0;
		
		task.materialConditions.byAreas.defaultMaterial = 
				std::make_shared<IsotropicMaterial>(rho0, lambda0, mu0);
		Task::MaterialCondition::ByAreas::Inhomogenity newMaterial;
		newMaterial.area = std::make_shared<AxisAlignedBoxArea>(
				Real3({-10, 0.5 - 1e-5, -10}), Real3({10, 10, 10}));
		newMaterial.material = std::make_shared<IsotropicMaterial>(rho, lambda, mu);
		task.materialConditions.byAreas.materials.push_back(newMaterial);
		
		task.cubicGrid.borderSize = 3;
		task.cubicGrid.h = {2.0/49, 1.0/99};
		task.cubicGrid.cubics = {{0, {{50, 100}, {0, 0}}}};
		
		task.globalSettings.numberOfSnaps = 0;
		task.globalSettings.requiredTime = 0.24;
		
		Task::InitialCondition::Wave wave;
		wave.waveType = Waves::T::P_FORWARD;
		wave.direction = 1; // along y
		wave.quantity = PhysicalQuantities::T::Vy;
		wave.quantityValue = -2;
		Real3 min({ -1, 0.015, -1});
		Real3 max({ 4, 0.455, 1});
		wave.area = std::make_shared<AxisAlignedBoxArea>(min, max);
		task.initialCondition.waves.push_back(wave);
		
		
		Engine<2> engine(task);
		auto init = Wrapper::getMesh(engine)->pdeVars({25, 25});
		ASSERT_NE(linal::zeros(init), init);
		engine.run();
		
//		int rightNodeIndex = (int) (task.cubicGrid.sizes(1) * 0.7);
		auto reflect = Wrapper::getMesh(engine)->pdeVars({25, 25});
//		auto transfer = Wrapper::getMesh(engine)->pdeVars(
//		{task.cubicGrid.sizes(0) / 2, rightNodeIndex});

		real E0 = mu0 * (3 * lambda0 + 2 * mu0) / (lambda0 + mu0); // Young's modulus
		real Z0 = sqrt(E0 * rho0);                                 // acoustic impedance
		real E = mu * (3 * lambda + 2 * mu) / (lambda + mu);       // Young's modulus
		real Z = sqrt(E * rho);                                    // acoustic impedance

		ASSERT_NEAR(reflect.sigma(1, 1) / init.sigma(1, 1),
		            (Z - Z0) / (Z + Z0),
		            1e-2);
		ASSERT_NEAR(reflect.velocity(1) / init.velocity(1),
		            (Z0 - Z) / (Z + Z0),
		            1e-2);

//		ASSERT_NEAR(transfer.sigma(1, 1) / init.sigma(1, 1),
//		            2 * Z / (Z + Z0),
//		            1e-2);
//		ASSERT_NEAR(transfer.velocity(1) / init.velocity(1),
//		            2 * Z0 / (Z + Z0),
//		            1e-2);


	}
}
