#include <lib/util/StringUtils.hpp>
#include <lib/Engine.hpp>
#include <lib/util/Area.hpp>

using namespace gcm;

inline Task parseTaskLayers() {
	Task task;

	task.modelId = Models::T::ELASTIC3D;
	task.materialId = Materials::T::ORTHOTROPIC;
	task.gridId = Grids::T::CUBIC;
	task.snapshottersId = {Snapshotters::T::VTK};

	task.cubicGrid.borderSize = 2;

	int debugDecrease = 3;
	real X = 0.016, Y = 0.016, Z = 0.004;
	task.cubicGrid.lengths = {X, Y, Z};
	task.cubicGrid.sizes = {151 / debugDecrease, 151 / debugDecrease, 101 / debugDecrease};

	Statement statement;
	OrthotropicMaterial defaultMaterial(1.6, 
			{163944.8, 3767.9, 3767.9,
					   8875.6, 2899.1,
							   8875.6,
									   4282.6, 4282.6, 4282.6});
	defaultMaterial.anglesOfRotation = {0, 0, M_PI / 4};
	statement.materialConditions.defaultMaterial = 
			std::make_shared<OrthotropicMaterial>(defaultMaterial);
	
	auto area90 = std::make_shared<AxisAlignedBoxArea>
			(Real3({-10, -10, Z / 2}), Real3({10, 10, 10}));
	auto rotated90 = defaultMaterial;
	rotated90.anglesOfRotation = {0, 0, -M_PI / 4};
	statement.materialConditions.materials.push_back({area90, 
			std::make_shared<OrthotropicMaterial>(rotated90)});

	
	statement.globalSettings.CourantNumber = 1.0;
	statement.globalSettings.numberOfSnaps = 200 / debugDecrease;
	statement.globalSettings.stepsPerSnap = 5;

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

	real fractureR = Z / 4;
//	Statement::Fracture fracture;
//	fracture.direction = 2;
//	fracture.index = task.cubicGrid.sizes.at(2) / 2;
//	fracture.area = std::make_shared<SphereArea>(Z / 2, Real3({X / 2, Y / 2, Z / 2}));
//	fracture.values = {
//		{PhysicalQuantities::T::Szz, [] (real) {return 0; }},
//		{PhysicalQuantities::T::Syz, [] (real) {return 0; }},
//		{PhysicalQuantities::T::Sxz, [] (real) {return 0; }}
//	};
//	statement.fractures.push_back(fracture);

	// quantities to snapshot
	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {PhysicalQuantities::T::PRESSURE};

	auto srcArea = std::make_shared<SphereArea>(fractureR, Real3({X/2, Y/2, Z}));
	real frequency = 10e+6, T = 1.0 / frequency;
	real A = -1e+6;
	borderCondition.area = srcArea;
	borderCondition.values = {
		{PhysicalQuantities::T::Szz, [A, T](real t) {return (t < T) ? A : 0; }},
		{PhysicalQuantities::T::Syz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Sxz, [] (real) {return 0; }}
	};
	statement.cubicGridBorderConditions.push_back(borderCondition);
	
	task.statements.push_back(statement);
	return task;
}


inline Task parseTaskCgalAnisotropy() {
	Task task;

	task.modelId = Models::T::ELASTIC3D;
	task.materialId = Materials::T::ISOTROPIC;
	task.gridId = Grids::T::CGAL;
	task.snapshottersId = {Snapshotters::T::VTK};

	task.cgal3DGrid.spatialStep = 0.002;
//	task.cgal3DGrid.detectSharpEdges = true;
	task.cgal3DGrid.polyhedronFileName = "meshes/layers_with_fracture.off";
	
	Statement statement;
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	statement.materialConditions.defaultMaterial =
	        std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);

	statement.globalSettings.CourantNumber = 0.5;

	statement.globalSettings.numberOfSnaps = 2;
	statement.globalSettings.stepsPerSnap = 1;

	Statement::BorderCondition borderConditionBottom;
	borderConditionBottom.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({0.01, 0.01, -10}), Real3({0.15, 0.15, 0.001}));
	borderConditionBottom.type = BorderConditions::T::FIXED_FORCE;
	borderConditionBottom.values = {
		[] (real) { return 0; },
		[] (real) { return 0; },
		[] (real) { return 0; }
	};
	
	Statement::BorderCondition borderConditionUp;
	borderConditionUp.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({0.01, 0.01, 0.039}), Real3({0.15, 0.15, 10}));
	borderConditionUp.type = BorderConditions::T::FIXED_FORCE;
	borderConditionUp.values = {
		[] (real) { return 0; },
		[] (real) { return 0; },
		[] (real) { return 0; }
	};
	
	Statement::BorderCondition sourceUp;
	sourceUp.area = std::make_shared<SphereArea>(0.01, Real3({0.08, 0.08, 0.04}));
	sourceUp.type = BorderConditions::T::FIXED_FORCE;
	sourceUp.values = {
		[] (real) { return 0; },
		[] (real) { return 0; },
		[] (real t) { return (t < 0.5) ? -1 : 0; }
	};
	
	Statement::BorderCondition fracture;
	fracture.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({0.02, 0.02, 0.01}), Real3({0.14, 0.14, 0.03}));
	fracture.type = BorderConditions::T::FIXED_FORCE;
	fracture.values = {
		[] (real) { return 0; },
		[] (real) { return 0; },
		[] (real) { return 0; }
	};

	statement.borderConditions = {borderConditionBottom,
	                              borderConditionUp,
	                              sourceUp};

	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {
		PhysicalQuantities::T::PRESSURE,
	};

	task.statements.push_back(statement);
	return task;
}
