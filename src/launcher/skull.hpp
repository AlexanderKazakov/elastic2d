#include <libgcm/util/task/Task.hpp>
#include <libgcm/util/math/Area.hpp>

using namespace gcm;


inline Task skullCommon() {
	Task task;
	task.calculationBasis = {
			1, 0, 0,
			0, 1, 0,
			0, 0, 1};
	
	task.globalSettings.dimensionality = 3;
	task.globalSettings.gridId = Grids::T::SIMPLEX;
	task.globalSettings.snapshottersId = { Snapshotters::T::VTK };
	task.globalSettings.CourantNumber = 1;
//	task.globalSettings.numberOfSnaps = 1000;
//	task.globalSettings.stepsPerSnap = 15;
	task.globalSettings.numberOfSnaps = 100;
	task.globalSettings.stepsPerSnap = 2;
	
	task.simplexGrid.mesher = Task::SimplexGrid::Mesher::INM_MESHER;
	task.simplexGrid.fileName = "meshes/coarse/skull_part.out";
//	task.simplexGrid.fileName = "meshes/coarse/skull-homogeneous.out";
//	task.simplexGrid.fileName = "meshes/coarse/mesh-aneurysm.out";
//	task.simplexGrid.fileName = "meshes/coarse/mesh-coarse.out";
//	task.simplexGrid.fileName = "meshes/refined/mesh-refined.out";
	task.simplexGrid.scale = 10;
	
	task.materialConditions.type = Task::MaterialCondition::Type::BY_BODIES;
	
	Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 1;
	pressure.area = std::make_shared<SphereArea>(1, Real3({0, 7.5, 142}));
	task.initialCondition.quantities.push_back(pressure);
	
	task.vtkSnapshotter.quantitiesToSnap = { PhysicalQuantities::T::PRESSURE };
	return task;
}


inline Task skullAcousticHomogeneous() {
	Task task = skullCommon();
	
	task.bodies = {
			{1, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {}}},
	};
	
	task.contactCondition.defaultCondition = ContactConditions::T::SLIDE;
	
	auto connectiveTissue = std::make_shared<IsotropicMaterial>(0.916,  1.886, 0, 0, 0, 1, 1.585);
	task.materialConditions.byBodies.bodyMaterialMap = {
			{1, connectiveTissue},
	};
	
	Task::BorderCondition freeBorder;
	freeBorder.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({-100, -100, 132}), Real3({100, 100, 1000}));
	freeBorder.type = BorderConditions::T::FIXED_FORCE;
	freeBorder.values = {
		[] (real) { return 0; },
	};
//	Task::BorderCondition source;
//	source.area = std::make_shared<SphereArea>(2, Real3({7, 3, 146.5}));
//	source.type = BorderConditions::T::FIXED_FORCE;
//	real tau = 0.5;
//	real omega = 2 * M_PI / tau;
//	source.values = {
//		[=] (real t) {
//			t -= 2 * tau;
//			return sin(omega * t) * exp(-t*t / ( 2 * tau*tau)); }
//	};
	task.borderConditions = {
		freeBorder,
//		source,
	};
	
	return task;
}


inline Task skullAcoustic() {
	Task task = skullCommon();
	
	task.bodies = {
			{1, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {}}},
			{2, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {}}},
			{3, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {}}},
			{4, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {}}},
			{5, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {}}},

//			{1, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {Odes::T::MAXWELL_VISCOSITY}}},
//			{2, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {Odes::T::MAXWELL_VISCOSITY}}},
//			{3, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {Odes::T::MAXWELL_VISCOSITY}}},
//			{4, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {}}},
//			{5, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {Odes::T::MAXWELL_VISCOSITY}}},
	};
	
	task.contactCondition.defaultCondition = ContactConditions::T::SLIDE;
	
	auto connectiveTissue = std::make_shared<IsotropicMaterial>(0.916,  1.886, 0, 0, 0, 1, 1.585);
	auto muscles          = std::make_shared<IsotropicMaterial>(1.041,  2.648, 0, 0, 0, 2, 0.878);
	auto cerebrum         = std::make_shared<IsotropicMaterial>(1.030,  2.475, 0, 0, 0, 3, 1.293);
	auto bones            = std::make_shared<IsotropicMaterial>(1.904,  7.854, 0, 0, 0, 4, 0);
	auto vessels          = std::make_shared<IsotropicMaterial>(1.066,  2.784, 0, 0, 0, 5, 1.288);
	task.materialConditions.byBodies.bodyMaterialMap = {
			{1, connectiveTissue},
			{2, muscles},
			{3, cerebrum},
			{4, bones},
			{5, vessels},
	};
	
	Task::BorderCondition freeBorder;
	freeBorder.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({-100, -100, 132}), Real3({100, 100, 1000}));
	freeBorder.type = BorderConditions::T::FIXED_FORCE;
	freeBorder.values = {
		[] (real) { return 0; },
	};
	Task::BorderCondition source;
	source.area = std::make_shared<SphereArea>(2, Real3({7, 3, 146.5}));
	source.type = BorderConditions::T::FIXED_FORCE;
	real tau = 0.5;
	real omega = 2 * M_PI / tau;
	source.values = {
		[=] (real t) {
			t -= 2 * tau;
			return sin(omega * t) * exp(-t*t / ( 2 * tau*tau)); }
	};
	task.borderConditions = {freeBorder, source};
	
	return task;
}


inline Task skull() {
	Task task = skullCommon();
	
	task.bodies = {
			{1, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}},
			{2, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}},
			{3, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}},
			{4, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}},
			{5, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}},
	};
	
	task.contactCondition.defaultCondition = ContactConditions::T::ADHESION;
	
	auto connectiveTissue = std::make_shared<IsotropicMaterial>(0.916,  1.415, 0.236, 0, 0, 1, 1.585);
	auto muscles          = std::make_shared<IsotropicMaterial>(1.041,  1.968, 0.331, 0, 0, 2, 0.878);
	auto cerebrum         = std::make_shared<IsotropicMaterial>(1.030,  1.856, 0.309, 0, 0, 3, 1.293);
	auto bones            = std::make_shared<IsotropicMaterial>(1.904,  5.891, 0.982, 0, 0, 4, 0);
	auto vessels          = std::make_shared<IsotropicMaterial>(1.066,  2.088, 0.348, 0, 0, 5, 1.288);
	task.materialConditions.byBodies.bodyMaterialMap = {
			{1, connectiveTissue},
			{2, muscles},
			{3, cerebrum},
			{4, bones},
			{5, vessels},
	};
	
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
	real tau = 0.5;
	real omega = 2 * M_PI / tau;
	source.values = {
		[] (real) { return 0; },
		[] (real) { return 0; },
		[=] (real t) {
			t -= 2 * tau;
			return sin(omega * t) * exp(-t*t / ( 2 * tau*tau)); }
	};
	task.borderConditions = {freeBorder, source};
	
	return task;
}
