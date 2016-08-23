#include <chrono>

#include <launcher/getopt_wrapper.hpp>
#include <launcher/skull.hpp>


using namespace gcm;

Task parseTaskCgal2d();
Task parseTaskCgal3d();
Task parseTask2d();
Task parseTask3d();
Task parseTaskCubicAcoustic();
Task parseTaskContact2D();


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
	else {
		LOG_FATAL("Invalid task file");
		return -1;
	}

	try {
		auto t1 = std::chrono::high_resolution_clock::now();
		Engine(task).run();
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
	
	task.bodies = {
		{0, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC}},
		{1, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC}}
	};
	
	task.simplexGrid.spatialStep = 0.2;
	
	Task::SimplexGrid::Body::Border body1border = {
		{3, 3}, {-3, 3}, {-3, -3}, {3, -3}
	};
	Task::SimplexGrid::Body::Border body2border = {
		{3, 3}, {9, 3}, {9, -3}, {3, -3}
	};
	task.simplexGrid.bodies = {
		Task::SimplexGrid::Body({0, body1border, {} }),
		Task::SimplexGrid::Body({1, body2border, {} })
	};
	
	Statement statement;
	statement.materialConditions.type = Statement::MaterialCondition::Type::BY_BODIES;
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	const auto material1 = std::make_shared<IsotropicMaterial>(rho, lambda, mu);
	const auto material2 = std::make_shared<IsotropicMaterial>(4 * rho, lambda, mu);
	statement.materialConditions.byBodies.bodyMaterialMap = {
		{0, material1},
		{1, material2}
	};
	
	statement.globalSettings.CourantNumber = 1;
	statement.globalSettings.numberOfSnaps = 50;
	statement.globalSettings.stepsPerSnap = 3;
	
	Statement::InitialCondition::Wave wave;
	wave.waveType = Waves::T::P_BACKWARD;
	wave.direction = 0;
	wave.quantity = PhysicalQuantities::T::PRESSURE;
	wave.quantityValue = 1;
	wave.area = std::make_shared<AxisAlignedBoxArea>(
				Real3({5, -10, -10}), Real3({7, 10, 10}));
	statement.initialCondition.waves.push_back(wave);
	
	Statement::BorderCondition borderConditionAll;
	borderConditionAll.area = std::make_shared<InfiniteArea>();
	borderConditionAll.type = BorderConditions::T::FIXED_FORCE;
	borderConditionAll.values = {
		[] (real) { return 0; },
//		[] (real) { return 0; }
	};
	
	Statement::BorderCondition borderConditionLeft;
	borderConditionLeft.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({-10, -10, -10}), Real3({-2.999, 10, 10}));
	borderConditionLeft.type = BorderConditions::T::FIXED_FORCE;
	borderConditionLeft.values = {
		[] (real) { return 0; },
		[] (real t) { return (t < 1) ? -1 : 0; }
	};
	
	Statement::BorderCondition borderConditionMid;
	borderConditionMid.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({-2.5, -2.5, -10}), Real3({0.5, 0.5, 10}));
	borderConditionMid.type = BorderConditions::T::FIXED_VELOCITY;
	borderConditionMid.values = {
		[] (real) { return 0; },
		[] (real) { return 0; }
	};
	
	statement.borderConditions = {borderConditionAll,
	                              /*borderConditionLeft,
	                              borderConditionMid*/};
	
	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {
		PhysicalQuantities::T::PRESSURE,
//		PhysicalQuantities::T::Sxx,
//		PhysicalQuantities::T::Sxy,
//		PhysicalQuantities::T::Syy
	};

	task.statements.push_back(statement);
	return task;
}


Task parseTaskCgal3d() {
	Task task;
	
	task.globalSettings.dimensionality = 3;
	task.globalSettings.gridId = Grids::T::SIMPLEX;
	task.globalSettings.snapshottersId = {Snapshotters::T::VTK};
	
	task.bodies = {{0, {Materials::T::ISOTROPIC, Models::T::ELASTIC}}};
	
	task.simplexGrid.mesher = Task::SimplexGrid::Mesher::CGAL_MESHER;
	task.simplexGrid.spatialStep = 0.2;
	task.simplexGrid.fileName = "meshes/icosahedron.off";
//	task.simplexGrid.detectSharpEdges = true;
//	task.simplexGrid.polyhedronFileName = "meshes/cube.off";
	
	Statement statement;
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	statement.materialConditions.byAreas.defaultMaterial =
	        std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);
	
	statement.globalSettings.CourantNumber = 1;
	
	statement.globalSettings.numberOfSnaps = 50;
	statement.globalSettings.stepsPerSnap = 1;
	
	Statement::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 0.5;
	pressure.area = std::make_shared<SphereArea>(0.5, Real3({0, 0, 0}));
//	pressure.area = std::make_shared<SphereArea>(0.2, Real3({0.5, 0.5, 0.5})); // for cube
	statement.initialCondition.quantities.push_back(pressure);

	Statement::BorderCondition borderConditionAll;
	borderConditionAll.area = std::make_shared<InfiniteArea>();
	borderConditionAll.type = BorderConditions::T::FIXED_FORCE;
	borderConditionAll.values = {
		[] (real) { return 0; },
		[] (real) { return 0; },
		[] (real) { return 0; }
	};
	
	Statement::BorderCondition borderConditionLeft;
	borderConditionLeft.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({-10, -10, -10}), Real3({0.01, 10, 10}));
	borderConditionLeft.type = BorderConditions::T::FIXED_FORCE;
	borderConditionLeft.values = {
		[] (real) { return 0; },
		[] (real) { return 0; },
		[] (real t) { return (t < 0.5) ? -1 : 0; }
	};
	
	Statement::BorderCondition borderConditionRight;
	borderConditionRight.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({0.99, -10, -10}), Real3({10, 10, 10}));
	borderConditionRight.type = BorderConditions::T::FIXED_VELOCITY;
	borderConditionRight.values = {
		[] (real) { return 0; },
		[] (real) { return 0; },
		[] (real) { return 0; }
	};
	
	statement.borderConditions = {borderConditionAll,
	                              /*borderConditionLeft,
	                              borderConditionRight*/};

	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {
		PhysicalQuantities::T::PRESSURE,
	};

	task.statements.push_back(statement);
	return task;
}


Task parseTaskCgal2d() {
	Task task;
	
	task.globalSettings.dimensionality = 2;
	task.globalSettings.gridId = Grids::T::SIMPLEX;
	task.globalSettings.snapshottersId = {Snapshotters::T::VTK};
	
	task.bodies = {
//			{0, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC}},
//			{1, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC}},
			{0, {Materials::T::ISOTROPIC, Models::T::ELASTIC}},
//			{1, {Materials::T::ISOTROPIC, Models::T::ELASTIC}}
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
	
	Statement statement;
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	statement.materialConditions.byAreas.defaultMaterial =
	        std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);
	
	statement.globalSettings.CourantNumber = 1;
	
	statement.globalSettings.numberOfSnaps = 50;
	statement.globalSettings.stepsPerSnap = 3;
	
	Statement::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 0.5;
	pressure.area = std::make_shared<SphereArea>(0.5, Real3({0, 6, 0}));
	statement.initialCondition.quantities.push_back(pressure);
	
//	Statement::InitialCondition::Wave wave;
//	wave.waveType = Waves::T::P_FORWARD;
//	wave.direction = 1;
//	wave.quantity = PhysicalQuantities::T::PRESSURE;
//	wave.quantityValue = 1;
//	wave.area = std::make_shared<AxisAlignedBoxArea>(Real3({-10, 0, -10}), Real3({10, 1, 10}));
//	statement.initialCondition.waves.push_back(wave);
	
	Statement::BorderCondition borderConditionAll;
	borderConditionAll.area = std::make_shared<InfiniteArea>();
	borderConditionAll.type = BorderConditions::T::FIXED_FORCE;
	borderConditionAll.values = {
		[] (real) { return 0; },
		[] (real) { return 0; }
	};
	
	Statement::BorderCondition borderConditionLeft;
	borderConditionLeft.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({-10, -10, -10}), Real3({-2.999, 10, 10}));
	borderConditionLeft.type = BorderConditions::T::FIXED_FORCE;
	borderConditionLeft.values = {
		[] (real) { return 0; },
		[] (real t) { return (t < 1) ? -1 : 0; }
	};
	
	Statement::BorderCondition borderConditionMid;
	borderConditionMid.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({-2.5, -2.5, -10}), Real3({0.5, 0.5, 10}));
	borderConditionMid.type = BorderConditions::T::FIXED_VELOCITY;
	borderConditionMid.values = {
		[] (real) { return 0; },
		[] (real) { return 0; }
	};
	
	statement.borderConditions = {borderConditionAll,
	                              borderConditionLeft,
	                              borderConditionMid};
	
	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {
		PhysicalQuantities::T::PRESSURE,
//		PhysicalQuantities::T::Sxx,
//		PhysicalQuantities::T::Sxy,
//		PhysicalQuantities::T::Syy
	};

	task.statements.push_back(statement);
	return task;
}


Task parseTask2d() {
	Task task;
	
	task.globalSettings.dimensionality = 2;
	task.globalSettings.gridId = Grids::T::CUBIC;
	task.globalSettings.snapshottersId = {Snapshotters::T::VTK};
	
	task.bodies = {{0, {Materials::T::ISOTROPIC, Models::T::ELASTIC}}};
	
	task.cubicGrid.borderSize = 2;
	task.cubicGrid.lengths = {4, 2, 1};
	task.cubicGrid.sizes = {100, 50, 1};

	Statement statement;
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	statement.materialConditions.byAreas.defaultMaterial =
	        std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);

	statement.globalSettings.CourantNumber = 0.9;

	statement.globalSettings.numberOfSnaps = 20;
	statement.globalSettings.stepsPerSnap = 1;

	Statement::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 10.0;
	pressure.area = std::make_shared<SphereArea>(0.2, Real3({2, 1, 0}));
	statement.initialCondition.quantities.push_back(pressure);

	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {PhysicalQuantities::T::PRESSURE};

	task.statements.push_back(statement);
	return task;
}


Task parseTask3d() {
	Task task;
	
	task.globalSettings.dimensionality = 3;
	task.globalSettings.gridId = Grids::T::CUBIC;
	task.globalSettings.snapshottersId = {Snapshotters::T::VTK};
	
	task.bodies = {{0, {Materials::T::ISOTROPIC, Models::T::ELASTIC}}};
	
	task.cubicGrid.borderSize = 2;
	task.cubicGrid.lengths = {4, 2, 1};
	task.cubicGrid.sizes = {200, 100, 50};

	Statement statement;
	real rho = 4;
//	statement.materialConditions.byAreas.defaultMaterial =
//	        std::make_shared<IsotropicMaterial>(rho, 2, 1, 1, 1);
	statement.materialConditions.byAreas.defaultMaterial = std::shared_ptr<OrthotropicMaterial>(
			new OrthotropicMaterial(rho, {360, 70, 70, 180, 70, 90, 10, 10, 10}, 1, 1));

	statement.globalSettings.CourantNumber = 1;

	statement.globalSettings.numberOfSnaps = 100;
	statement.globalSettings.stepsPerSnap = 1;

	Statement::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 10.0;
	pressure.area = std::make_shared<SphereArea>(0.2, Real3({2, 1, 0.5}));
	statement.initialCondition.quantities.push_back(pressure);

	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {
			PhysicalQuantities::T::PRESSURE,
			PhysicalQuantities::T::Sxx,
			PhysicalQuantities::T::Sxy,
			PhysicalQuantities::T::Sxz,
			PhysicalQuantities::T::Syy,
			PhysicalQuantities::T::Syz,
			PhysicalQuantities::T::Szz
	};

	task.statements.push_back(statement);
	return task;
}


inline Task parseTaskCubicAcoustic() {
	Task task;
	
	task.globalSettings.dimensionality = 2;
	task.globalSettings.gridId = Grids::T::CUBIC;
	task.globalSettings.snapshottersId = {Snapshotters::T::VTK};
	
	task.bodies = {{0, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC}}};
	
	task.cubicGrid.borderSize = 2;
	task.cubicGrid.lengths = {2, 1};//, 1};
	task.cubicGrid.sizes = {100, 50};//, 50};

	Statement statement;
	statement.materialConditions.byAreas.defaultMaterial =
		std::make_shared<IsotropicMaterial>(1, 1, 0);

	statement.globalSettings.CourantNumber = 1;

	statement.globalSettings.numberOfSnaps = 80;
	statement.globalSettings.stepsPerSnap = 1;

	Statement::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 10.0;
	pressure.area = std::make_shared<SphereArea>(0.2, Real3({1, 0.5, 0}));
	statement.initialCondition.quantities.push_back(pressure);

	
	Statement::CubicGridBorderCondition borderCondition;
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
			(Real3({-10, -10, -10}), Real3({10, 1e-5, 10}));
	borderCondition.values = {
			{PhysicalQuantities::T::PRESSURE, [] (real) {return 0; }},
	};
	statement.cubicGridBorderConditions.push_back(borderCondition);
	
	
	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {
			PhysicalQuantities::T::PRESSURE,
	};

	task.statements.push_back(statement);
	return task;
}
