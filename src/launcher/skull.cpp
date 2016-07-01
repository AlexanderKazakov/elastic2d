#include <lib/util/StringUtils.hpp>
#include <lib/Engine.hpp>
#include <lib/util/Area.hpp>

using namespace gcm;


inline Task skull() {
	Task task;

	task.modelId = Models::T::ELASTIC3D;
	task.materialId = Materials::T::ISOTROPIC;
	task.gridId = Grids::T::CGAL;
	task.snapshottersId = {
			Snapshotters::T::VTK
	};

	
	task.cgal3DGrid.mesher = Task::Cgal3DGrid::Mesher::INM_MESHER;
	task.cgal3DGrid.spatialStep = 0.5;
	task.cgal3DGrid.fileName = "meshes/coarse/mesh-coarse.out";
	task.cgal3DGrid.scale = 10;
	
	
	Statement statement;
	statement.materialConditions.type = Statement::MaterialCondition::Type::BY_CELLS;
	auto connectiveTissue = std::make_shared<IsotropicMaterial>(1.008, 2.187, 0.0911, 0, 0, 1);
	auto muscles          = std::make_shared<IsotropicMaterial>(1.041, 1.765, 0.4413, 0, 0, 2);
	auto cerebrum         = std::make_shared<IsotropicMaterial>(1.030, 2.284, 0.0952, 0, 0, 3);
	auto bones            = std::make_shared<IsotropicMaterial>(1.672, 5.197, 2.6773, 0, 0, 4);
	auto vessels          = std::make_shared<IsotropicMaterial>(1.063, 2.517, 0.1049, 0, 0, 5);
	statement.materialConditions.materialMap = {
			{1, {connectiveTissue, 100} },
			{2, {muscles,          1} },
			{3, {cerebrum,         5} },
			{4, {bones,            10} },
			{5, {vessels,          1000} },
	};
	
	statement.globalSettings.CourantNumber = 1;
	statement.globalSettings.numberOfSnaps = 100;
	statement.globalSettings.stepsPerSnap = 5;

	
	Statement::BorderCondition freeBorder;
	freeBorder.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({-100, -100, 132}), Real3({100, 100, 1000}));
	freeBorder.type = BorderConditions::T::FIXED_FORCE;
	freeBorder.values = {
		[] (real) { return 0; },
		[] (real) { return 0; },
		[] (real) { return 0; }
	};
	
	
//	Statement::BorderCondition source;
//	source.area = std::make_shared<SphereArea>(2, Real3({-7, 3, 146.5}));
//	source.type = BorderConditions::T::FIXED_FORCE;
//	source.values = {
//		[] (real) { return 0; },
//		[] (real) { return 0; },
//		[] (real /*t*/) { return -1; }
////		[] (real t) { return sin(omega * t) * exp(-t*t / ( 2 * tau)); }
//	};
	
	statement.borderConditions = {/*freeBorder, /*source*/};

	Statement::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 1;
	pressure.area = std::make_shared<SphereArea>(2, Real3({0, 0, 148}));
	statement.initialCondition.quantities.push_back(pressure);
	
	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {
		PhysicalQuantities::T::PRESSURE,
	};

	task.statements.push_back(statement);
	return task;
}
