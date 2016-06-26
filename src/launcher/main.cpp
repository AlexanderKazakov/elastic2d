#include <chrono>

#include <launcher/getopt_wrapper.hpp>
#include <lib/util/StringUtils.hpp>
#include <lib/Engine.hpp>
#include <lib/util/Area.hpp>

#include <launcher/diploma_tasks.cpp>

const int NUMBER_OF_SENSOR_POSITIONS_ALONG_AXIS = 10;

using namespace gcm;

Task parseTaskCgal2d();
Task parseTaskCgal3d();
Task parseTask2d();
Task parseTask3d();
Task parseTaskSeismo();
Task parseTaskCagi2d();
Task parseTaskCagi3d();


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	USE_AND_INIT_LOGGER("gcm.main");
	
	std::string taskId;
	getopt_wrapper(argc, argv, taskId);
	
	Task task;
	if      (taskId == "cgal2d"    ) { task = parseTaskCgal2d(); }
	else if (taskId == "cgal3d"    ) { task = parseTaskCgal3d(); }
	else if (taskId == "seismo"    ) { task = parseTaskSeismo(); }
	else if (taskId == "cubic"     ) { task = parseTask3d();     }
	else if (taskId == "inverse"   ) { task = parse3D(); }
	else if (taskId == "layers"    ) { task = parseTaskLayers(); }
	else if (taskId == "cgalani"   ) { task = parseTaskCgalAnisotropy(); }
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


Task parseTaskCgal3d() {
	Task task;

	task.modelId = Models::T::ELASTIC3D;
	task.materialId = Materials::T::ISOTROPIC;
	task.gridId = Grids::T::CGAL;
	task.snapshottersId = {Snapshotters::T::VTK};

	task.cgal3DGrid.spatialStep = 0.1;
	task.cgal3DGrid.polyhedronFileName = "meshes/icosahedron.off";
//	task.cgal3DGrid.detectSharpEdges = true;
//	task.cgal3DGrid.polyhedronFileName = "meshes/cube.off";
	
	Statement statement;
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	statement.materialConditions.defaultMaterial =
	        std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);

	statement.globalSettings.CourantNumber = 0.5;

	statement.globalSettings.numberOfSnaps = 70;
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

	task.modelId = Models::T::ELASTIC2D;
	task.materialId = Materials::T::ISOTROPIC;
	task.gridId = Grids::T::CGAL;
	task.snapshottersId = {Snapshotters::T::VTK};

	task.cgal2DGrid.spatialStep = 0.2;
	task.cgal2DGrid.movable = false;
	
	Task::Cgal2DGrid::Body::Border outer = {
		{3, 3}, {-3, 3}, {-3, -3}, {3, -3}, {2, 1},
	};
	task.cgal2DGrid.bodies = {
		Task::Cgal2DGrid::Body(outer,
				{/*{{-2, -1}, {-1, 0}, {0, -1}, {-1, -2}}, {{1, 1}, {1, 2}, {2, 2}, {2, 1}}*/}),
//		Task::Cgal2DGrid::Body({{-2, 5}, {2, 5}, {0, 7}}, {})
	};

	Statement statement;
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	statement.materialConditions.defaultMaterial =
	        std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);

	statement.globalSettings.CourantNumber = 1;

	statement.globalSettings.numberOfSnaps = 40;
	statement.globalSettings.stepsPerSnap = 1;

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
		PhysicalQuantities::T::Sxx,
		PhysicalQuantities::T::Sxy,
		PhysicalQuantities::T::Syy
	};

	task.statements.push_back(statement);
	return task;
}


Task parseTask2d() {
	Task task;

	task.modelId = Models::T::ELASTIC2D;
	task.materialId = Materials::T::ISOTROPIC;
	task.gridId = Grids::T::CUBIC;
	task.snapshottersId = {Snapshotters::T::VTK};

	task.cubicGrid.borderSize = 2;
	task.cubicGrid.lengths = {4, 2, 1};
	task.cubicGrid.sizes = {100, 50, 1};

	Statement statement;
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	statement.materialConditions.defaultMaterial =
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

	task.modelId = Models::T::ELASTIC3D;
	task.materialId = Materials::T::ORTHOTROPIC;
	task.gridId = Grids::T::CUBIC;
	task.snapshottersId = {Snapshotters::T::VTK};

	task.cubicGrid.borderSize = 2;
	task.cubicGrid.lengths = {4, 2, 1};
	task.cubicGrid.sizes = {200, 100, 50};

	Statement statement;
	real rho = 4;
//	statement.materialConditions.defaultMaterial =
//	        std::make_shared<IsotropicMaterial>(rho, 2, 1, 1, 1);
	statement.materialConditions.defaultMaterial = std::shared_ptr<OrthotropicMaterial>(
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


Task parseTaskSeismo() {
	Task task;

	task.modelId = Models::T::ELASTIC2D;
	task.materialId = Materials::T::ISOTROPIC;
	task.gridId = Grids::T::CUBIC;
	task.snapshottersId = {
		Snapshotters::T::BIN2DSEISM,
		Snapshotters::T::VTK
	};

	task.globalSettings.forceSequence = true;

	task.cubicGrid.borderSize = 3;
	task.cubicGrid.lengths = {1, 1};
	task.cubicGrid.sizes = {50, 50};

	Statement statement;
	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.binary2DSeismograph.quantityToWrite = PhysicalQuantities::T::PRESSURE;
	statement.globalSettings.CourantNumber = 1.0; // number from Courant–Friedrichs–Lewy
	                                              // condition

	real rho = 4;                                 // default density
	real lambda = 2;                              // default Lame parameter
	real mu = 1;                                  // default Lame parameter
	statement.materialConditions.defaultMaterial =
	        std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);

	statement.globalSettings.numberOfSnaps = 50;
	statement.globalSettings.stepsPerSnap = 1;

	Statement::InitialCondition::Wave wave;
	wave.waveType = Waves::T::P_FORWARD;
	wave.direction = 1;
	wave.quantity = PhysicalQuantities::T::PRESSURE;
	wave.quantityValue = MPI::COMM_WORLD.Get_rank() + 1;
	wave.area = std::make_shared<AxisAlignedBoxArea>(Real3({-1, 0.2, -1}), Real3({3, 0.5, 3}));
	statement.initialCondition.waves.push_back(wave);

	statement.id = "0000";
	task.statements.push_back(statement);

	Statement::CubicGridBorderCondition borderCondition;
	// y right free border
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
			(Real3({-10, 0.99, -10}), Real3({10, 10, 10}));
	borderCondition.values = {
		{PhysicalQuantities::T::Sxy, [] (real) {return 0; }},
		{PhysicalQuantities::T::Syy, [] (real) {return 0; }}
	};
	statement.cubicGridBorderConditions.push_back(borderCondition);

	statement.id = "0001";
	task.statements.push_back(statement);

	return task;
}


Task parseTaskCagi2d() {
	Task task;

	task.modelId = Models::T::ELASTIC2D;
	task.materialId = Materials::T::ISOTROPIC;
	task.gridId = Grids::T::CUBIC;
	task.snapshottersId = {
			Snapshotters::T::VTK, 
			Snapshotters::T::DETECTOR
	};

	task.cubicGrid.borderSize = 2;
	task.globalSettings.forceSequence = true;

	int debugDecrease = 3;

	real X = 0.016, Y = 0.004;
	real sensorSize = 0.003;
	real sourceSize = 0.003;
	task.cubicGrid.lengths = {X, Y};
	task.cubicGrid.sizes = {151 / debugDecrease, 101 / debugDecrease};

	Statement statement;
	real rho = 1e+3;
	real lambda = 3e+10;
	real mu = 2e+10;
	statement.materialConditions.defaultMaterial =
	        std::make_shared<IsotropicMaterial>(rho, lambda, mu);
	statement.globalSettings.CourantNumber = 1.0;
	statement.globalSettings.numberOfSnaps = 251 / debugDecrease;
	statement.globalSettings.stepsPerSnap = 1;

	Statement::CubicGridBorderCondition borderCondition;
	// y bottom free border
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
			(Real3({-10, -10, -10}), Real3({10, 1e-5, 10}));
	borderCondition.values = {
		{PhysicalQuantities::T::Sxy, [] (real) {return 0; }},
		{PhysicalQuantities::T::Syy, [] (real) {return 0; }}
	};
	statement.cubicGridBorderConditions.push_back(borderCondition);
	// y up free border
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
			(Real3({-10, Y - 1e-5, -10}), Real3({10, 10, 10}));
	borderCondition.values = {
		{PhysicalQuantities::T::Sxy, [] (real) {return 0; }},
		{PhysicalQuantities::T::Syy, [] (real) {return 0; }}
	};
	statement.cubicGridBorderConditions.push_back(borderCondition);

	Statement::Fracture fracture;
	fracture.direction = 1;
	fracture.index = task.cubicGrid.sizes.at(1) / 2;
	fracture.area = std::make_shared<SphereArea>(Y / 2, Real3({X / 2, Y / 2, 0}));
	fracture.values = {
		{PhysicalQuantities::T::Sxy, [] (real) {return 0; }},
		{PhysicalQuantities::T::Syy, [] (real) {return 0; }}
	};
	statement.fractures.push_back(fracture);

	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {PhysicalQuantities::T::PRESSURE};
	statement.detector.quantities = {PhysicalQuantities::T::Vy};

	real spaceToMoveSensor = X - (sensorSize + sourceSize);
	real shiftAtTime = spaceToMoveSensor / NUMBER_OF_SENSOR_POSITIONS_ALONG_AXIS;
	std::cout << "shiftAtTime = " << shiftAtTime << std::endl;
	int counter = 0;
	for (int i = 0; i < NUMBER_OF_SENSOR_POSITIONS_ALONG_AXIS; i++) {
		real initSensorPosition = i * shiftAtTime;
		real endSensorPosition = initSensorPosition + sensorSize;
		real initSrcPosition = endSensorPosition;
		real endSrcPosition = initSrcPosition + sourceSize;
		// y up sensor
		auto srcArea = std::make_shared<AxisAlignedBoxArea>
				(Real3({initSrcPosition, Y - 1e-5, -10}), Real3({endSrcPosition, 10, 10}));
		real frequency = 10e+6, T = 1.0 / frequency;
		real A = -1e+6;
		borderCondition.area = srcArea;
		borderCondition.values = {
			{PhysicalQuantities::T::Syy, [A, T](real t) {return (t < T) ? A : 0; }},
			{PhysicalQuantities::T::Sxy, [] (real) {return 0; }}
		};
		if (statement.cubicGridBorderConditions.size() == 3) {statement.cubicGridBorderConditions.pop_back(); }
		statement.cubicGridBorderConditions.push_back(borderCondition);
		auto sensorArea = std::make_shared<AxisAlignedBoxArea>
				(Real3({initSensorPosition, Y - 1e-5, -10}), Real3({endSensorPosition, 10, 10}));
		statement.detector.area = sensorArea;

		statement.id = StringUtils::toString(i, 4);
		if (counter % MPI::COMM_WORLD.Get_size() == MPI::COMM_WORLD.Get_rank()) {
//			if (i == 3) {
				task.statements.push_back(statement);
//			}
		}
		counter++;
	}
	return task;
}


Task parseTaskCagi3d() {
	Task task;

	task.modelId = Models::T::ELASTIC3D;
	task.materialId = Materials::T::ISOTROPIC;
	task.gridId = Grids::T::CUBIC;
	task.snapshottersId = {/*Snapshotters::T::VTK,*/ Snapshotters::T::DETECTOR};

	task.cubicGrid.borderSize = 2;
	task.globalSettings.forceSequence = true;

	int debugDecrease = 3;
	real X = 0.016, Y = 0.016, Z = 0.004;
	real sensorSizeX = 0.003;
	real sourceSizeX = 0.003;
	real commonSizeY = 0.006;
	task.cubicGrid.lengths = {X, Y, Z};
	task.cubicGrid.sizes = {151 / debugDecrease, 151 / debugDecrease, 101 / debugDecrease};

	Statement statement;
	real rho = 1e+3;
	real lambda = 3e+10;
	real mu = 2e+10;
	statement.materialConditions.defaultMaterial =
	        std::make_shared<IsotropicMaterial>(rho, lambda, mu);
	statement.globalSettings.CourantNumber = 1.0;
	statement.globalSettings.numberOfSnaps = 251 / debugDecrease;
	statement.globalSettings.stepsPerSnap = 1;

	Statement::CubicGridBorderCondition borderCondition;
	// z bottom free border
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
			(Real3({-10, -10, -10}), Real3({10, 10, 1e-5}));
	borderCondition.values = {
		{PhysicalQuantities::T::Szz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Syz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Sxz, [] (real) {return 0; }}
	};
	statement.cubicGridBorderConditions.push_back(borderCondition);
	// z up free border
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
			(Real3({-10, -10, Z - 1e-5}), Real3({10, 10, 10}));
	borderCondition.values = {
		{PhysicalQuantities::T::Szz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Syz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Sxz, [] (real) {return 0; }}
	};
	statement.cubicGridBorderConditions.push_back(borderCondition);

	Statement::Fracture fracture;
	fracture.direction = 2;
	fracture.index = task.cubicGrid.sizes.at(2) / 2;
	fracture.area = std::make_shared<SphereArea>(Z / 2, Real3({X / 2, Y / 2, Z / 2}));
	fracture.values = {
		{PhysicalQuantities::T::Szz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Syz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Sxz, [] (real) {return 0; }}
	};
	statement.fractures.push_back(fracture);

	// quantities to snapshot
	statement.detector.quantities = {PhysicalQuantities::T::Vz};
	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {PhysicalQuantities::T::PRESSURE};

	real spaceToMoveSensorX = X - (sensorSizeX + sourceSizeX);
	real shiftAtTimeX = spaceToMoveSensorX / NUMBER_OF_SENSOR_POSITIONS_ALONG_AXIS;
	std::cout << "shiftAtTimeX = " << shiftAtTimeX << std::endl;
	real spaceToMoveSensorY = Y - commonSizeY;
	real shiftAtTimeY = spaceToMoveSensorY / NUMBER_OF_SENSOR_POSITIONS_ALONG_AXIS;
	std::cout << "shiftAtTimeY = " << shiftAtTimeY << std::endl;
	int counter = 0;
	for (int i = 0; i < NUMBER_OF_SENSOR_POSITIONS_ALONG_AXIS; i++) {
		for (int j = 0; j < NUMBER_OF_SENSOR_POSITIONS_ALONG_AXIS; j++) {
			// detector x
			real initSensorPositionX = i * shiftAtTimeX;
			real endSensorPositionX = initSensorPositionX + sensorSizeX;
			real initSrcPositionX = endSensorPositionX;
			real endSrcPositionX = initSrcPositionX + sourceSizeX;
			// detector y
			real initSensorPositionY = j * shiftAtTimeY;
			real endSensorPositionY = initSensorPositionY + commonSizeY;
			real initSrcPositionY = initSensorPositionY;
			real endSrcPositionY = endSensorPositionY;
			// z up sensor
			auto srcArea = std::make_shared<AxisAlignedBoxArea>
					(Real3({initSrcPositionX, initSrcPositionY, Z - 1e-5}),
					 Real3({endSrcPositionX, endSrcPositionY, 10}));
			real frequency = 10e+6, T = 1.0 / frequency;
			real A = -1e+6;
			borderCondition.area = srcArea;
			borderCondition.values = {
				{PhysicalQuantities::T::Szz, [A, T](real t) {return (t < T) ? A : 0; }},
				{PhysicalQuantities::T::Syz, [] (real) {return 0; }},
				{PhysicalQuantities::T::Sxz, [] (real) {return 0; }}
			};
			if (statement.cubicGridBorderConditions.size() == 3) {
				statement.cubicGridBorderConditions.pop_back();
			}
			statement.cubicGridBorderConditions.push_back(borderCondition);
			auto sensorArea = std::make_shared<AxisAlignedBoxArea>
					(Real3({initSensorPositionX, initSensorPositionY, Z - 1e-5}),
					 Real3({endSensorPositionX, endSensorPositionY, 10}));
			statement.detector.area = sensorArea;

			statement.id = StringUtils::toString(j, 2) + StringUtils::toString(i, 2);
			if (counter % MPI::COMM_WORLD.Get_size() == MPI::COMM_WORLD.Get_rank()) {
				if (i == 5 /*&& j == 5*/) {
					task.statements.push_back(statement);
				}
			}
			counter++;
		}
	}
	return task;
}
