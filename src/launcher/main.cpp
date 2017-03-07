#include <chrono>

#include <libgcm/engine/EngineFactory.hpp>

#include <launcher/getopt_wrapper.hpp>
#include <launcher/skull.hpp>
#include <launcher/ndi.hpp>
#include <launcher/tor.hpp>


using namespace gcm;

Task parseTaskCgal2d();
Task parseTaskCgal3d();
Task parseTask2d();
Task parseTask3d();
Task parseTaskCubicAcoustic();
Task parseTaskContact2D();
Task parseTaskTmp();
Task parseTaskCube();

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	USE_AND_INIT_LOGGER("gcm.main");
	
	std::string taskId;
	getopt_wrapper(argc, argv, taskId);
	
	Task task;
	if      (taskId == "cgal2d"    ) { task = parseTaskCgal2d(); }
	else if (taskId == "cgal3d"    ) { task = parseTaskCgal3d(); }
	else if (taskId == "cubic2d"   ) { task = parseTask2d();     }
	else if (taskId == "cubic3d"   ) { task = parseTask3d();     }
	else if (taskId == "acoustic"  ) { task = parseTaskCubicAcoustic(); }
	else if (taskId == "skull"     ) { task = skull(); }
	else if (taskId == "skullAcs"  ) { task = skullAcoustic(); }
	else if (taskId == "contact"   ) { task = parseTaskContact2D(); }
	else if (taskId == "tmp"       ) { task = parseTaskTmp(); }
	else if (taskId == "cube"      ) { task = parseTaskCube(); }
	else if (taskId == "ndi"       ) { task = ndi(); }
	else if (taskId == "tor"       ) { task = torAcoustic(); }
	else {
		LOG_FATAL("Invalid task file");
		return -1;
	}
	
	try {
		auto t1 = std::chrono::high_resolution_clock::now();
		createEngine(task)->run();
		auto t2 = std::chrono::high_resolution_clock::now();
		SUPPRESS_WUNUSED(t1); SUPPRESS_WUNUSED(t2);

		LOG_INFO("Time of calculation, microseconds = " << 
			std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count());
	} catch (Exception e) {
		LOG_FATAL(e.what());
	}

	MPI_Finalize();
	return 0;
}


Task parseTaskContact2D() {
	Task task;
	
	task.globalSettings.dimensionality = 2;
	task.globalSettings.gridId = Grids::T::SIMPLEX;
	task.globalSettings.snapshottersId = {Snapshotters::T::VTK};
	task.contactCondition.defaultCondition = ContactConditions::T::SLIDE;
	task.globalSettings.CourantNumber = 1;
	task.globalSettings.numberOfSnaps = 150;
	task.globalSettings.stepsPerSnap = 3;
	
	task.bodies = {
		{0, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {/*Odes::T::MAXWELL_VISCOSITY*/}}},
		{1, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {/*Odes::T::MAXWELL_VISCOSITY*/}}}
	};
	
	task.simplexGrid.spatialStep = 0.2;
	
	Task::SimplexGrid::Body::Border body1border = {
		{3, 3}, {-3, 3}, {-3, -3}, {3, -3}
	};
	Task::SimplexGrid::Body::Border body2border = {
		{3, 3}, {9, 3}, {9, -3}, {3, -3}
	};
//	Task::SimplexGrid::Body::Border body2border = {
//		{-4, 4}, {0, 6}, {4, 4}
//	};
	task.simplexGrid.bodies = {
		Task::SimplexGrid::Body({0, body1border, {} }),
		Task::SimplexGrid::Body({1, body2border, {} })
	};
	
	task.materialConditions.type = Task::MaterialCondition::Type::BY_BODIES;
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	const auto material1 = std::make_shared<IsotropicMaterial>(4 * rho, lambda, mu, 0, 0, 0, 0);
	const auto material2 = std::make_shared<IsotropicMaterial>(rho, lambda, mu, 0, 0, 0, 1e+9);
	task.materialConditions.byBodies.bodyMaterialMap = {
		{0, material1},
		{1, material2}
	};
	
	Task::BorderCondition borderConditionAll;
	borderConditionAll.area = std::make_shared<InfiniteArea>();
	borderConditionAll.type = BorderConditions::T::FIXED_FORCE;
	borderConditionAll.values = {
		[] (real) { return 0; },
//		[] (real) { return 0; }
	};
	
	Task::BorderCondition borderConditionLeft;
	borderConditionLeft.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({-10, -10, -10}), Real3({-2.999, 10, 10}));
	borderConditionLeft.type = BorderConditions::T::FIXED_FORCE;
	borderConditionLeft.values = {
		[] (real) { return 0; },
		[] (real t) { return (t < 1) ? -1 : 0; }
	};
	
	Task::BorderCondition borderConditionMid;
	borderConditionMid.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({-2.5, -2.5, -10}), Real3({0.5, 0.5, 10}));
	borderConditionMid.type = BorderConditions::T::FIXED_VELOCITY;
	borderConditionMid.values = {
		[] (real) { return 0; },
		[] (real) { return 0; }
	};
	
	task.borderConditions = {borderConditionAll,
	                              /*borderConditionLeft,
	                              borderConditionMid*/};
	
	Task::InitialCondition::Wave wave;
	wave.waveType = Waves::T::P_BACKWARD;
	wave.direction = 0;
	wave.quantity = PhysicalQuantities::T::PRESSURE;
	wave.quantityValue = 1;
	wave.area = std::make_shared<AxisAlignedBoxArea>(
				Real3({5, -10, -10}), Real3({7, 10, 10}));
	task.initialCondition.waves.push_back(wave);
	
	task.vtkSnapshotter.quantitiesToSnap = {
		PhysicalQuantities::T::PRESSURE,
//		PhysicalQuantities::T::Sxx,
//		PhysicalQuantities::T::Sxy,
//		PhysicalQuantities::T::Syy
	};
	
	return task;
}


Task parseTaskCgal3d() {
	Task task;
	task.calculationBasis = {
			1, 0, 0,
			0, 1, 0,
			0, 0, 1};
	
	task.globalSettings.dimensionality = 3;
	task.globalSettings.gridId = Grids::T::SIMPLEX;
	task.globalSettings.snapshottersId = {Snapshotters::T::VTK};
	
	task.bodies = {{0, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {/*Odes::T::MAXWELL_VISCOSITY*/}}}};
	
	task.simplexGrid.mesher = Task::SimplexGrid::Mesher::CGAL_MESHER;
	task.simplexGrid.spatialStep = 0.2;
	task.simplexGrid.fileName = "meshes/icosahedron.off";
//	task.simplexGrid.spatialStep = 0.1;
//	task.simplexGrid.detectSharpEdges = true;
//	task.simplexGrid.fileName = "meshes/cube.off";
	
	real rho = 4;
	real lambda = 1; // 2;
	real mu = 1;
	task.materialConditions.type = Task::MaterialCondition::Type::BY_BODIES;
	const auto material = std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1, 1, 1);
	task.materialConditions.byBodies.bodyMaterialMap = {
		{0, material},
	};
	
	task.globalSettings.CourantNumber = 1;
	
	task.globalSettings.numberOfSnaps = 50;
	task.globalSettings.stepsPerSnap = 3;
	
	Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 0.5;
	pressure.area = std::make_shared<SphereArea>(0.5, Real3({0, 0, 0}));
//	pressure.area = std::make_shared<SphereArea>(0.2, Real3({0.5, 0.5, 0.5})); // for cube
	task.initialCondition.quantities.push_back(pressure);

	Task::BorderCondition borderConditionAll;
	borderConditionAll.area = std::make_shared<InfiniteArea>();
	borderConditionAll.type = BorderConditions::T::FIXED_FORCE;
	borderConditionAll.values = {
		[] (real) { return 0; },
		[] (real) { return 0; },
		[] (real) { return 0; }
	};
	
	task.borderConditions = {borderConditionAll,
	                              /*borderConditionLeft,
	                              borderConditionRight*/};
	task.vtkSnapshotter.quantitiesToSnap = {
		PhysicalQuantities::T::PRESSURE,
	};
	return task;
}


Task parseTaskCgal2d() {
	Task task;
	real phi = M_PI / 4;
	task.calculationBasis = {
			cos(phi), -sin(phi),
			sin(phi),  cos(phi),
	};
	
	task.globalSettings.dimensionality = 2;
	task.globalSettings.gridId = Grids::T::SIMPLEX;
	task.globalSettings.snapshottersId = {Snapshotters::T::VTK};
	
	task.bodies = {
//			{0, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {}}},
//			{1, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {}}},
			{0, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}},
//			{1, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}}
	};
	
	task.simplexGrid.spatialStep = 0.2;
	
	Task::SimplexGrid::Body::Border outer = {
		{3, 3}, {-3, 3}, {-3, -3}, {3, -3}, {2, 2},
	};
	std::vector<Task::SimplexGrid::Body::Border> inner = {{
		{-2, -1}, {-1, 0}, {0, -1}, {-1, -2}
	}};
	task.simplexGrid.bodies = {
		Task::SimplexGrid::Body({0, outer, inner}),
//		Task::SimplexGrid::Body({1, { {-2, 5}, {2, 5}, {0, 7} }, {}})
	};
	
//	task.contactCondition.defaultCondition = ContactConditions::T::SLIDE;
	task.contactCondition.defaultCondition = ContactConditions::T::ADHESION;
	
	real rho = 4;
	real lambda = 1; // 2;
	real mu = 1;
	task.materialConditions.type = Task::MaterialCondition::Type::BY_BODIES;
	const auto material = std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);
	task.materialConditions.byBodies.bodyMaterialMap = {
		{0, material},
	};
	
	task.globalSettings.CourantNumber = 1;
	task.globalSettings.numberOfSnaps = 150;
	task.globalSettings.stepsPerSnap = 1;
	
	Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 0.5;
	pressure.area = std::make_shared<SphereArea>(0.5, Real3({0, 6, 0}));
	task.initialCondition.quantities.push_back(pressure);
	
//	Task::InitialCondition::Wave wave;
//	wave.waveType = Waves::T::P_FORWARD;
//	wave.direction = 1;
//	wave.quantity = PhysicalQuantities::T::PRESSURE;
//	wave.quantityValue = 1;
//	wave.area = std::make_shared<AxisAlignedBoxArea>(Real3({-10, 0, -10}), Real3({10, 1, 10}));
//	task.initialCondition.waves.push_back(wave);
	
	Task::BorderCondition borderConditionAll;
	borderConditionAll.area = std::make_shared<InfiniteArea>();
	borderConditionAll.type = BorderConditions::T::FIXED_FORCE;
	borderConditionAll.values = {
		[] (real) { return 0; },
		[] (real) { return 0; }
	};
	
	Task::BorderCondition borderConditionLeft;
	borderConditionLeft.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({-10, -10, -10}), Real3({-2.999, 10, 10}));
	borderConditionLeft.type = BorderConditions::T::FIXED_FORCE;
	borderConditionLeft.values = {
		[] (real) { return 0; },
		[] (real t) { return (t < 1) ? -1 : 0; }
	};
	
	Task::BorderCondition borderConditionMid;
	borderConditionMid.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({-2.5, -2.5, -10}), Real3({0.5, 0.5, 10}));
	borderConditionMid.type = BorderConditions::T::FIXED_FORCE;
	borderConditionMid.values = {
		[] (real) { return 0; },
		[] (real) { return 0; }
	};
	
	task.borderConditions = {borderConditionAll,
	                              borderConditionLeft,
	                              borderConditionMid};
	
	
	task.vtkSnapshotter.quantitiesToSnap = {
		PhysicalQuantities::T::PRESSURE,
		PhysicalQuantities::T::Sxx,
		PhysicalQuantities::T::Sxy,
		PhysicalQuantities::T::Syy
	};
	
	return task;
}


Task parseTask2d() {
	Task task;
	
	task.globalSettings.dimensionality = 2;
	task.globalSettings.gridId = Grids::T::CUBIC;
	task.globalSettings.snapshottersId = {Snapshotters::T::VTK};
	
	task.bodies = {{0, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}}};
	
	task.cubicGrid.borderSize = 2;
	task.cubicGrid.h = {4.0/99, 2.0/49};
	task.cubicGrid.cubics = {{0, {{100, 50}, {0, 0}}}};
	
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	task.materialConditions.byAreas.defaultMaterial =
	        std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);
	
	task.globalSettings.CourantNumber = 0.9;
	
	task.globalSettings.numberOfSnaps = 20;
	task.globalSettings.stepsPerSnap = 1;
	
	Task::InitialCondition::Wave wave;
	wave.waveType = Waves::T::P_FORWARD;
	wave.direction = 0;
	wave.quantity = PhysicalQuantities::T::PRESSURE;
	wave.quantityValue = 1;
	wave.area = std::make_shared<AxisAlignedBoxArea>(
				Real3({1, -10, -10}), Real3({2, 10, 10}));
	task.initialCondition.waves.push_back(wave);
	
//	Task::InitialCondition::Quantity pressure;
//	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
//	pressure.value = 10.0;
//	pressure.area = std::make_shared<SphereArea>(0.2, Real3({2, 1, 0}));
//	task.initialCondition.quantities.push_back(pressure);
	
	task.vtkSnapshotter.quantitiesToSnap = {PhysicalQuantities::T::PRESSURE};
	return task;
}


Task parseTask3d() {
	Task task;
	
	task.globalSettings.dimensionality = 3;
	task.globalSettings.gridId = Grids::T::CUBIC;
	task.globalSettings.snapshottersId = {Snapshotters::T::VTK};
	
	task.bodies = {{0, {Materials::T::ORTHOTROPIC, Models::T::ELASTIC, {}}}};
	
	task.cubicGrid.borderSize = 2;
	task.cubicGrid.h = {4.0/99, 2.0/49, 1.0/19};
	task.cubicGrid.cubics = {{0, {{100, 50, 20}, {0, 0, 0}}}};

	real rho = 4;
//	task.materialConditions.byAreas.defaultMaterial =
//	        std::make_shared<IsotropicMaterial>(rho, 2, 1, 1, 1);
	task.materialConditions.byAreas.defaultMaterial = std::shared_ptr<OrthotropicMaterial>(
			new OrthotropicMaterial(rho, {360, 70, 70, 180, 70, 90, 10, 10, 10}, 1, 1));

	task.globalSettings.CourantNumber = 1;

	task.globalSettings.numberOfSnaps = 100;
	task.globalSettings.stepsPerSnap = 1;

	Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 10.0;
	pressure.area = std::make_shared<SphereArea>(0.2, Real3({2, 1, 0.5}));
	task.initialCondition.quantities.push_back(pressure);

	
	task.vtkSnapshotter.quantitiesToSnap = {
			PhysicalQuantities::T::PRESSURE,
			PhysicalQuantities::T::Sxx,
			PhysicalQuantities::T::Sxy,
			PhysicalQuantities::T::Sxz,
			PhysicalQuantities::T::Syy,
			PhysicalQuantities::T::Syz,
			PhysicalQuantities::T::Szz
	};

	
	return task;
}


inline Task parseTaskCubicAcoustic() {
	Task task;
	
	task.globalSettings.dimensionality = 2;
	task.globalSettings.gridId = Grids::T::CUBIC;
	task.globalSettings.snapshottersId = {Snapshotters::T::VTK};
	
	task.bodies = {{0, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {}}}};
	
	task.cubicGrid.borderSize = 2;
	task.cubicGrid.h = {1, 1};
	task.cubicGrid.cubics = {{0, {{100, 50}, {0, 0}}}};
	
	task.materialConditions.byAreas.defaultMaterial =
		std::make_shared<IsotropicMaterial>(1, 1, 0);
	
	task.globalSettings.CourantNumber = 1;
	
	task.globalSettings.numberOfSnaps = 80;
	task.globalSettings.stepsPerSnap = 1;
	
	Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 10.0;
	pressure.area = std::make_shared<SphereArea>(10, Real3({50, 25, 0}));
	task.initialCondition.quantities.push_back(pressure);
	
	
	auto area = std::make_shared<AxisAlignedBoxArea>
			(Real3({-1e+3, -1e+3, -1e+3}), Real3({1e+3, 1e-5, 1e+3}));
	Task::CubicBorderCondition::Values freeBorder = {
		{PhysicalQuantities::T::PRESSURE, [] (real) {return 0; }},
	};
	task.cubicBorderConditions = {
		{0, {
			 {1, area, freeBorder},
		 }},
	};
	
	
	task.vtkSnapshotter.quantitiesToSnap = {
			PhysicalQuantities::T::PRESSURE,
	};
	
	return task;
}


inline Task parseTaskTmp() {
	Task task;
	
	task.globalSettings.dimensionality = 2;
	task.globalSettings.gridId = Grids::T::SIMPLEX;
	task.globalSettings.snapshottersId = {Snapshotters::T::VTK};
	task.contactCondition.defaultCondition = ContactConditions::T::SLIDE;
	task.globalSettings.CourantNumber = 1.0;
	task.globalSettings.requiredTime = 14;
//	task.globalSettings.numberOfSnaps = 100;
	task.globalSettings.stepsPerSnap = 3;
	
//	real phi = M_PI / 4;
//	task.calculationBasis = {
//			cos(phi), -sin(phi),
//			sin(phi),  cos(phi),
//	};
	task.calculationBasis = {
			1, 0,
			0, 1
	};
	
	task.bodies = {
		{0, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {}}},
//		{1, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {}}}
	};
	
	task.simplexGrid.spatialStep = 0.1;
	
	Task::SimplexGrid::Body::Border body1border = {
		{0, 3}, {-3, 3}, {-3, -3}, {3, -3}
	};
//	Task::SimplexGrid::Body::Border body2border = {
//		{3, 3}, {9, 3}, {9, -3}, {3, -3}
//	};
	task.simplexGrid.bodies = {
		Task::SimplexGrid::Body({0, body1border, {} }),
//		Task::SimplexGrid::Body({1, body2border, {} })
	};
	
	task.materialConditions.type = Task::MaterialCondition::Type::BY_BODIES;
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	const auto material1 = std::make_shared<IsotropicMaterial>(rho, lambda, mu, 0, 0, 0, 0);
	const auto material2 = std::make_shared<IsotropicMaterial>(rho, lambda, mu, 0, 0, 0, 1e+9);
	task.materialConditions.byBodies.bodyMaterialMap = {
		{0, material1},
		{1, material2}
	};
	
	Task::BorderCondition borderConditionAll;
	borderConditionAll.area = std::make_shared<InfiniteArea>();
	borderConditionAll.type = BorderConditions::T::FIXED_FORCE;
	borderConditionAll.values = {
		[] (real) { return 0; },
	};
	Task::BorderCondition borderConditionLeft;
	borderConditionLeft.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({-10, -2.99999, -10}), Real3({-2.99999, 2.99999, 10}));
	borderConditionLeft.type = BorderConditions::T::FIXED_FORCE;
	borderConditionLeft.values = {
		[] (real t) { return (t < 1) ? 1 : 0; }
	};
	task.borderConditions = {borderConditionAll,
	                         borderConditionLeft};
	
	
	task.vtkSnapshotter.quantitiesToSnap = {
		PhysicalQuantities::T::PRESSURE,
	};
	return task;
}



Task parseTaskCube() {
	Task task;
	task.calculationBasis = {
			1, 0, 0,
			0, 1, 0,
			0, 0, 1};
//	linal::createLocalBasis(linal::normalize(Real3({1, 1, 1})));
//	task.calculationBasis = 
	
	task.globalSettings.dimensionality = 3;
	task.globalSettings.gridId = Grids::T::SIMPLEX;
	task.globalSettings.snapshottersId = {Snapshotters::T::VTK};
	task.bodies = {{0, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}}};
//	task.bodies = {{0, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {}}}};
	
	const int N = 2;
	task.simplexGrid.mesher = Task::SimplexGrid::Mesher::CGAL_MESHER;
	task.simplexGrid.spatialStep = 0.025 * N;
	task.simplexGrid.detectSharpEdges = true;
	task.simplexGrid.fileName = "meshes/cube.off";
	
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	task.materialConditions.type = Task::MaterialCondition::Type::BY_BODIES;
	const auto material = std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);
	task.materialConditions.byBodies.bodyMaterialMap = {
		{0, material},
	};
	
	task.globalSettings.CourantNumber = 1;
	task.globalSettings.numberOfSnaps = 200 / N;
	task.globalSettings.stepsPerSnap = 1;
	
//	Task::InitialCondition::Quantity pressure;
//	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
//	pressure.value = 0.5;
//	pressure.area = std::make_shared<SphereArea>(0.2, Real3({0.5, 0.5, 0.5}));
//	task.initialCondition.quantities.push_back(pressure);

	Task::BorderCondition borderConditionAll;
	borderConditionAll.area = std::make_shared<InfiniteArea>();
	borderConditionAll.type = BorderConditions::T::FIXED_FORCE;
	borderConditionAll.values = {
		[] (real) { return 0; },
		[] (real) { return 0; },
		[] (real) { return 0; }
	};
	
	Task::BorderCondition borderConditionLeft;
	borderConditionLeft.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({-10, -10, -10}), Real3({0.01, 10, 10}));
	borderConditionLeft.type = BorderConditions::T::FIXED_FORCE;
	borderConditionLeft.values = {
		[] (real) { return 0; },
		[] (real) { return 0; },
		[] (real t) { return (t < 0.5) ? -1 : 0; }
	};
	
	Task::BorderCondition borderConditionRight;
	borderConditionRight.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({0.99, -10, -10}), Real3({10, 10, 10}));
	borderConditionRight.type = BorderConditions::T::FIXED_VELOCITY;
	borderConditionRight.values = {
		[] (real) { return 0; },
		[] (real) { return 0; },
		[] (real) { return 0; }
	};
	
	task.borderConditions = {borderConditionAll,
	                         borderConditionLeft,
	                       /*borderConditionRight*/};
	
	task.vtkSnapshotter.quantitiesToSnap = {
		PhysicalQuantities::T::PRESSURE,
	};
	return task;
}

