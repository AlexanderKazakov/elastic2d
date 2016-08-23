#include <lib/Engine.hpp>
#include <lib/util/Area.hpp>

using namespace gcm;


inline Task skullAcoustic() {
	Task task;
	
	task.globalSettings.dimensionality = 3;
	task.globalSettings.gridId = Grids::T::SIMPLEX;
	task.globalSettings.snapshottersId = { Snapshotters::T::VTK };
	
	task.bodies = {
			{1, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC}},
			{2, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC}},
			{3, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC}},
			{4, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC}},
			{5, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC}},
	};
	
	
	task.simplexGrid.mesher = Task::SimplexGrid::Mesher::INM_MESHER;
	task.simplexGrid.fileName = "meshes/refined/mesh-refined.out";
	task.simplexGrid.scale = 10;
	
	task.contactCondition.defaultCondition = ContactConditions::T::SLIDE;
	
	task.materialConditions.type = Task::MaterialCondition::Type::BY_BODIES;
	auto connectiveTissue = std::make_shared<IsotropicMaterial>(1.008,  2.369, 0, 0, 0, 1);
	auto muscles          = std::make_shared<IsotropicMaterial>(1.041,  2.648, 0, 0, 0, 2);
	auto cerebrum         = std::make_shared<IsotropicMaterial>(1.030,  2.475, 0, 0, 0, 3);
	auto bones            = std::make_shared<IsotropicMaterial>(1.672, 10.552, 0, 0, 0, 4);
	auto vessels          = std::make_shared<IsotropicMaterial>(1.063,  2.726, 0, 0, 0, 5);
	task.materialConditions.byBodies.bodyMaterialMap = {
			{1, connectiveTissue},
			{2, muscles},
			{3, cerebrum},
			{4, bones},
			{5, vessels},
	};
	
	
	task.globalSettings.CourantNumber = 1;
	task.globalSettings.numberOfSnaps = 1000;
	task.globalSettings.stepsPerSnap = 5;
	
	
	Task::BorderCondition freeBorder;
	freeBorder.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({-100, -100, 132}), Real3({100, 100, 1000}));
	freeBorder.type = BorderConditions::T::FIXED_FORCE;
	freeBorder.values = {
		[] (real) { return 0; },
	};
	
	
	Task::BorderCondition source;
	source.area = std::make_shared<SphereArea>(2, Real3({-7, 3, 146.5}));
	source.type = BorderConditions::T::FIXED_FORCE;
	real tau = 0.5;
	real omega = 2 * M_PI / tau;
	source.values = {
//		[=] (real t) { return (t < tau) ? 1 : 0; }
		[=] (real t) {
			t -= 2 * tau;
			return sin(omega * t) * exp(-t*t / ( 2 * tau*tau)); }
	};
	
	task.borderConditions = {freeBorder, source};
	
//	Task::InitialCondition::Quantity pressure;
//	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
//	pressure.value = 1;
//	pressure.area = std::make_shared<SphereArea>(2, Real3({0, 5, 147}));
//	task.initialCondition.quantities.push_back(pressure);
	
	
	task.vtkSnapshotter.quantitiesToSnap = {
		PhysicalQuantities::T::PRESSURE,
	};
	
	return task;
}


inline Task skull() {
	Task task;
	
	task.globalSettings.dimensionality = 3;
	task.globalSettings.gridId = Grids::T::SIMPLEX;
	task.globalSettings.snapshottersId = { Snapshotters::T::VTK };
	
	task.bodies = {
			{1, {Materials::T::ISOTROPIC, Models::T::ELASTIC}},
			{2, {Materials::T::ISOTROPIC, Models::T::ELASTIC}},
			{3, {Materials::T::ISOTROPIC, Models::T::ELASTIC}},
			{4, {Materials::T::ISOTROPIC, Models::T::ELASTIC}},
			{5, {Materials::T::ISOTROPIC, Models::T::ELASTIC}},
	};
	
	
	task.simplexGrid.mesher = Task::SimplexGrid::Mesher::INM_MESHER;
	task.simplexGrid.fileName = "meshes/coarse/mesh-coarse.out";
	task.simplexGrid.scale = 10;
	
	task.contactCondition.defaultCondition = ContactConditions::T::ADHESION;
	
	
//	task.materialConditions.type = Task::MaterialCondition::Type::BY_AREAS;
//	task.materialConditions.byAreas.defaultMaterial = 
//			std::make_shared<IsotropicMaterial>(1, 2, 1);
	
	task.materialConditions.type = Task::MaterialCondition::Type::BY_BODIES;
	auto connectiveTissue = std::make_shared<IsotropicMaterial>(1.008, 2.187, 0.0911, 0, 0, 1);
	auto muscles          = std::make_shared<IsotropicMaterial>(1.041, 1.765, 0.4413, 0, 0, 2);
	auto cerebrum         = std::make_shared<IsotropicMaterial>(1.030, 2.284, 0.0952, 0, 0, 3);
	auto bones            = std::make_shared<IsotropicMaterial>(1.672, 5.197, 2.6773, 0, 0, 4);
	auto vessels          = std::make_shared<IsotropicMaterial>(1.063, 2.517, 0.1049, 0, 0, 5);
	task.materialConditions.byBodies.bodyMaterialMap = {
			{1, connectiveTissue},
			{2, muscles},
			{3, cerebrum},
			{4, bones},
			{5, vessels},
	};
	
	task.globalSettings.CourantNumber = 1;
	task.globalSettings.numberOfSnaps = 1000;
	task.globalSettings.stepsPerSnap = 1;
	
	
	Task::BorderCondition freeBorder;
	freeBorder.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({-100, -100, 132}), Real3({100, 100, 1000}));
	freeBorder.type = BorderConditions::T::FIXED_FORCE;
	freeBorder.values = {
		[] (real) { return 0; },
		[] (real) { return 0; },
		[] (real) { return 0; }
	};
	
	
	Task::BorderCondition source;
	source.area = std::make_shared<SphereArea>(2, Real3({-7, 3, 146.5}));
	source.type = BorderConditions::T::FIXED_FORCE;
	source.values = {
		[] (real) { return 0; },
		[] (real) { return 0; },
		[] (real t) { return (t < 0.5) ? -1 : 0; }
//		[] (real t) { return sin(omega * t) * exp(-t*t / ( 2 * tau)); }
	};
	
	task.borderConditions = {freeBorder, source};
	
	Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 1;
	pressure.area = std::make_shared<SphereArea>(2, Real3({0, 5, 147}));
	task.initialCondition.quantities.push_back(pressure);
	
	
	task.vtkSnapshotter.quantitiesToSnap = {
		PhysicalQuantities::T::PRESSURE,
	};
	
	return task;
}
