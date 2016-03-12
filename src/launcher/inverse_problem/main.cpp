#include <lib/Engine.hpp>
#include <lib/rheology/models/Model.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/util/snapshot/VtkSnapshotter.hpp>
#include <lib/util/snapshot/Detector.hpp>
#include <lib/util/areas/areas.hpp>


using namespace gcm;

Task parseTaskCagi2d();
Task parseTaskCagi3d();

const int NUMBER_OF_SENSOR_POSITIONS_ALONG_AXIS = 10;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	USE_AND_INIT_LOGGER("gcm.main");
	try {
		Engine engine;
		engine.setSolver(new DefaultSolver<DefaultMesh<Elastic2DModel, CubicGrid>>());
		engine.addSnapshotter(new VtkSnapshotter<DefaultMesh<Elastic2DModel, CubicGrid>>());
		engine.addSnapshotter(new Detector<DefaultMesh<Elastic2DModel, CubicGrid>>());
		engine.initialize(parseTaskCagi2d());
		engine.run();
	} catch (Exception e) {
		LOG_FATAL(e.what());
	}

	MPI_Finalize();
	return 0;
}


Task parseTaskCagi2d() {
	Task task;
	task.dimensionality = 2;
	task.borderSize = 2;
	task.forceSequence = true;
	task.enableSnapshotting = true;

	real X = 0.016, Y = 0.004;
	real sensorSize = 0.003;
	real sourceSize = 0.003;
	task.lengthes = {X, Y, 1};
	task.sizes = {151, 101, 1};

	Statement statement;
	real rho = 1e+3;
	real lambda = 3e+10;
	real mu = 2e+10;
	statement.isotropicMaterial = IsotropicMaterial(rho, lambda, mu);
	statement.CourantNumber = 1.0;
	statement.numberOfSnaps = 251;
	statement.stepsPerSnap = 1;

	Statement::BorderCondition borderCondition;	
	// y bottom free border
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
		(linal::Vector3({-10, -10, -10}), linal::Vector3({10, 1e-5, 10}));
	borderCondition.values = {
		{PhysicalQuantities::T::Sxy, [](real){return 0;}},
		{PhysicalQuantities::T::Syy, [](real){return 0;}}
	};
	statement.borderConditions.push_back(borderCondition);
	// y up free border
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
		(linal::Vector3({-10, Y-1e-5, -10}), linal::Vector3({10, 10, 10}));
	borderCondition.values = {
		{PhysicalQuantities::T::Sxy, [](real){return 0;}},
		{PhysicalQuantities::T::Syy, [](real){return 0;}}
	};
	statement.borderConditions.push_back(borderCondition);
	
	Statement::Fracture fracture;
	fracture.direction = 1;
	fracture.coordinate = Y/2;
	fracture.area = std::make_shared<SphereArea>(Y/2, linal::Vector3({X/2, Y/2, 0}));
	fracture.values = {
		{PhysicalQuantities::T::Sxy, [](real){return 0;}},
		{PhysicalQuantities::T::Syy, [](real){return 0;}}
	};
	statement.fractures.push_back(fracture); 
	
	// quantities to snapshot
	statement.quantitiesToVtk = {};
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
			(linal::Vector3({initSrcPosition, Y - 1e-5, -10}),
			 linal::Vector3({endSrcPosition, 10, 10}));
		real frequency = 10e+6, T = 1.0 / frequency;
		real A = - 1e+6;
		borderCondition.area = srcArea;
		borderCondition.values = {
			{PhysicalQuantities::T::Syy, [A, T](real t){return (t < T) ? A : 0;}},
			{PhysicalQuantities::T::Sxy, [](real){return 0;}}
		};
		if (statement.borderConditions.size() == 3) statement.borderConditions.pop_back();
		statement.borderConditions.push_back(borderCondition);
		auto sensorArea = std::make_shared<AxisAlignedBoxArea>
			(linal::Vector3({initSensorPosition, Y - 1e-5, -10}),
			 linal::Vector3({endSensorPosition, 10, 10}));
		statement.detector.area = sensorArea;

		statement.id = StringUtils::toString(i, 4);
		if (counter % MPI::COMM_WORLD.Get_size() == MPI::COMM_WORLD.Get_rank()) {
			if (i == 5) {
				task.statements.push_back(statement);
			}
		}
		counter++;
	}
	return task;
}


Task parseTaskCagi3d() {
	Task task;
	task.dimensionality = 3;
	task.borderSize = 2;
	task.forceSequence = true;
	task.enableSnapshotting = true;

	real X = 0.016, Y = 0.016, Z = 0.004;
	real sensorSizeX = 0.003;
	real sourceSizeX = 0.003;
	real commonSizeY = 0.006;
	task.lengthes = {X, Y, Z};
	task.sizes = {151/2, 151/2, 101/2};

	Statement statement;
	real rho = 1e+3;
	real lambda = 3e+10;
	real mu = 2e+10;
	statement.isotropicMaterial = IsotropicMaterial(rho, lambda, mu);
	statement.CourantNumber = 1.0;
	statement.numberOfSnaps = 251/2;
	statement.stepsPerSnap = 1;

	Statement::BorderCondition borderCondition;	
	// z bottom free border
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
		(linal::Vector3({-10, -10, -10}), linal::Vector3({10, 10, 1e-5}));
	borderCondition.values = {
		{PhysicalQuantities::T::Szz, [](real){return 0;}},
		{PhysicalQuantities::T::Syz, [](real){return 0;}},
		{PhysicalQuantities::T::Sxz, [](real){return 0;}}
	};
	statement.borderConditions.push_back(borderCondition);
	// z up free border
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
		(linal::Vector3({-10, -10, Z-1e-5}), linal::Vector3({10, 10, 10}));
	borderCondition.values = {
		{PhysicalQuantities::T::Szz, [](real){return 0;}},
		{PhysicalQuantities::T::Syz, [](real){return 0;}},
		{PhysicalQuantities::T::Sxz, [](real){return 0;}}
	};
	statement.borderConditions.push_back(borderCondition);
	
	Statement::Fracture fracture;
	fracture.direction = 2;
	fracture.coordinate = Z/2;
	fracture.area = std::make_shared<SphereArea>(Z/2, linal::Vector3({X/2, Y/2, Z/2}));
	fracture.values = {
		{PhysicalQuantities::T::Szz, [](real){return 0;}},
		{PhysicalQuantities::T::Syz, [](real){return 0;}},
		{PhysicalQuantities::T::Sxz, [](real){return 0;}}
	};
	statement.fractures.push_back(fracture); 
	
	// quantities to snapshot
	statement.quantitiesToVtk = {};
	statement.detector.quantities = {PhysicalQuantities::T::Vz};
	
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
				(linal::Vector3({initSrcPositionX, initSrcPositionY, Z - 1e-5}),
				 linal::Vector3({endSrcPositionX, endSrcPositionY, 10}));
			real frequency = 10e+6, T = 1.0 / frequency;
			real A = - 1e+6;
			borderCondition.area = srcArea;
			borderCondition.values = {
				{PhysicalQuantities::T::Szz, [A, T](real t){return (t < T) ? A : 0;}},
				{PhysicalQuantities::T::Syz, [](real){return 0;}},
				{PhysicalQuantities::T::Sxz, [](real){return 0;}}
			};
			if (statement.borderConditions.size() == 3) statement.borderConditions.pop_back();
			statement.borderConditions.push_back(borderCondition);
			auto sensorArea = std::make_shared<AxisAlignedBoxArea>
				(linal::Vector3({initSensorPositionX, initSensorPositionY, Z - 1e-5}),
				 linal::Vector3({endSensorPositionX, endSensorPositionY, 10}));
			statement.detector.area = sensorArea;

			statement.id = StringUtils::toString(j, 2) + StringUtils::toString(i, 2);
			if (counter % MPI::COMM_WORLD.Get_size() == MPI::COMM_WORLD.Get_rank()) {
				if (i == 5 && j == 5)
					task.statements.push_back(statement);
			}
			counter++;
		}
	}
	return task;
}
