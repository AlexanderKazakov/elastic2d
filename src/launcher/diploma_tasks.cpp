#include <lib/util/StringUtils.hpp>
#include <lib/Engine.hpp>
#include <lib/util/Area.hpp>

using namespace gcm;

const real sizeX = 0.016, sizeY = 0.016, sizeZ = 0.004;
const int nLayers = 11;
const real layerWidth = sizeZ / nLayers;

inline std::shared_ptr<OrthotropicMaterial>
createMaterial(const real phi) {
	OrthotropicMaterial material(1.6, 
			{163944.8, 3767.9, 3767.9,
					   8875.6, 2899.1,
							   8875.6,
									   4282.6, 4282.6, 4282.6});
	material.anglesOfRotation = {0, 0, phi};
	return std::make_shared<OrthotropicMaterial>(material);
}


inline std::shared_ptr<AxisAlignedBoxArea>
createArea(const int layerNumber /* from 0 at the bottom */) {
	const real Zmin = layerWidth * layerNumber;
	const real Zmax = layerWidth * (layerNumber + 1);
	const real Xmin = -10, Ymin = -10;
	const real Xmax = 100 * sizeX;
	const real Ymax = 100 * sizeY;
	return std::make_shared<AxisAlignedBoxArea>
			(Real3({Xmin, Ymin, Zmin}), Real3({Xmax, Ymax, Zmax}));
}



inline Task parseTaskTriangles() {
	Task task;

	task.modelId = Models::T::ELASTIC2D;
	task.materialId = Materials::T::ISOTROPIC;
	task.gridId = Grids::T::CGAL;
	task.snapshottersId = {Snapshotters::T::VTK};

	task.cgal2DGrid.spatialStep = 0.002;
	task.cgal2DGrid.movable = false;
	
	Task::Cgal2DGrid::Body::Border outer = {
		{0, 0}, {0.16, 0}, {0.16, 0.04}, {0, 0.04}
	};
	task.cgal2DGrid.bodies = {
		Task::Cgal2DGrid::Body(outer, 
				{{ {0.06, 0.02}, {0.10, 0.02}, {0.08, 0.022} }}),
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

	Statement::BorderCondition borderConditionAll;
	borderConditionAll.area = std::make_shared<InfiniteArea>();
	borderConditionAll.type = BorderConditions::T::FIXED_FORCE;
	borderConditionAll.values = {
		[] (real) { return 0; },
		[] (real) { return 0; }
	};

	Statement::BorderCondition sourceUp;
	sourceUp.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({0.06, 0.03, -10}), Real3({0.10, 0.05, 10}));
	sourceUp.type = BorderConditions::T::FIXED_FORCE;
	sourceUp.values = {
		[] (real) { return 0; },
		[] (real t) { return (t < 1) ? -1 : 0; }
	};
	
	Statement::BorderCondition borderConditionMid;
	borderConditionMid.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({0.01, 0.01, -10}), Real3({0.15, 0.03, 10}));
	borderConditionMid.type = BorderConditions::T::FIXED_FORCE;
	borderConditionMid.values = {
		[] (real) { return 0; },
		[] (real) { return 0; }
	};
	
	statement.borderConditions = {borderConditionAll,
	                              sourceUp,
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




inline Task parseTaskLayers() {
	Task task;

	task.modelId = Models::T::ELASTIC3D;
	task.materialId = Materials::T::ORTHOTROPIC;
	task.gridId = Grids::T::CUBIC;
	task.snapshottersId = {
			Snapshotters::T::VTK,
			Snapshotters::T::SLICESNAP,
	};

	task.cubicGrid.borderSize = 2;

	int debugDecrease = 2;
	task.cubicGrid.lengths = {sizeX, sizeY, sizeZ};
	task.cubicGrid.sizes = {301 / debugDecrease, 301 / debugDecrease, 201 / debugDecrease};

	Statement statement;
	
	statement.materialConditions.defaultMaterial = createMaterial(0);
	statement.materialConditions.materials = {
			{createArea(0),  createMaterial( M_PI / 4)},
			{createArea(2),  createMaterial(-M_PI / 4)},
			{createArea(5),  createMaterial( M_PI / 2)},
			{createArea(8),  createMaterial(-M_PI / 4)},
			{createArea(10), createMaterial( M_PI / 4)},
			/*
			{createArea(11), createMaterial( M_PI / 4)},
			{createArea(13), createMaterial(-M_PI / 4)},
			{createArea(16), createMaterial( M_PI / 2)},
			{createArea(19), createMaterial(-M_PI / 4)},
			{createArea(21), createMaterial( M_PI / 4)},
			*/
	};
	
	statement.globalSettings.CourantNumber = 1.0;
	statement.globalSettings.numberOfSnaps = 2000 / debugDecrease;
	statement.globalSettings.stepsPerSnap = 15;

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
			(Real3({-10, -10, sizeZ - 1e-5}), Real3({10, 10, 10}));
	borderCondition.values = {
		{PhysicalQuantities::T::Szz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Syz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Sxz, [] (real) {return 0; }}
	};
	statement.cubicGridBorderConditions.push_back(borderCondition);

	real fractureR = sizeZ / 4;
	
	Statement::Fracture fracture;
	fracture.direction = 2;
	fracture.index = task.cubicGrid.sizes.at(2) / nLayers;
	fracture.area = std::make_shared<SphereArea>(fractureR, Real3({sizeX / 2, sizeY / 2, layerWidth}));
	fracture.values = {
		{PhysicalQuantities::T::Szz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Syz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Sxz, [] (real) {return 0; }}
	};
//	statement.fractures.push_back(fracture);
	
	
	// quantities to snapshot
	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {
			PhysicalQuantities::T::Szz,
			//PhysicalQuantities::T::Sxz,
			//PhysicalQuantities::T::Syz,
			PhysicalQuantities::T::PRESSURE,
	};
	
	auto srcArea = std::make_shared<SphereArea>(2*fractureR, Real3({sizeX/2, sizeY/2, sizeZ}));
	statement.detector.quantities = {PhysicalQuantities::T::Vz};
	statement.detector.area = srcArea;
	
	real frequency = 10e+6, T = 1.0 / frequency;
	real A = -1e+6;
	borderCondition.area = srcArea;
	borderCondition.values = {
		{PhysicalQuantities::T::Szz, [A, T](real t) {return (t < T) ? A : 0; }},
		{PhysicalQuantities::T::Syz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Sxz, [] (real) {return 0; }},
	};
	statement.cubicGridBorderConditions.push_back(borderCondition);
	
	task.statements.push_back(statement);
	return task;
}


inline Task parse2D() {
	Task task;

	task.modelId = Models::T::ELASTIC2D;
	task.materialId = Materials::T::ISOTROPIC;
	task.gridId = Grids::T::CUBIC;
	task.snapshottersId = {
			Snapshotters::T::VTK,
	};

	task.cubicGrid.borderSize = 2;

	int debugDecrease = 1;

	real X = 0.016, Y = 0.004;
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

	// y up sensor
	auto srcArea = std::make_shared<SphereArea>(Y / 2, Real3({X / 2, Y, 0}));
	real frequency = 10e+6, T = 1.0 / frequency;
	real A = -1e+6;
	borderCondition.area = srcArea;
	borderCondition.values = {
			{PhysicalQuantities::T::Syy, [A, T](real t) {return (t < T) ? A : 0; }},
			{PhysicalQuantities::T::Sxy, [] (real) {return 0; }}
	};
	
	statement.cubicGridBorderConditions.push_back(borderCondition);
	
	task.statements.push_back(statement);
	return task;
}


inline Task parse3D() {
	Task task;

	task.modelId = Models::T::ELASTIC3D;
	task.materialId = Materials::T::ISOTROPIC;
	task.gridId = Grids::T::CUBIC;
	task.snapshottersId = {
			Snapshotters::T::VTK,
	};

	task.cubicGrid.borderSize = 2;

	int debugDecrease = 2;

	real X = 0.016, Y = 0.016, Z = 0.004;
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
		{PhysicalQuantities::T::Sxz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Syz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Szz, [] (real) {return 0; }},
	};
	statement.cubicGridBorderConditions.push_back(borderCondition);
	// z up free border
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
			(Real3({-10, -10, Z - 1e-5}), Real3({10, 10, 10}));
	borderCondition.values = {
		{PhysicalQuantities::T::Sxz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Syz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Szz, [] (real) {return 0; }},	};
	statement.cubicGridBorderConditions.push_back(borderCondition);

	Statement::Fracture fracture;
	fracture.direction = 2;
	fracture.index = task.cubicGrid.sizes.at(2) / 2;
	fracture.area = std::make_shared<SphereArea>(Z / 2, Real3({X / 2, Y / 2, Z / 2}));
	fracture.values = {
		{PhysicalQuantities::T::Sxz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Syz, [] (real) {return 0; }},
		{PhysicalQuantities::T::Szz, [] (real) {return 0; }},
	};
	statement.fractures.push_back(fracture);

	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {PhysicalQuantities::T::PRESSURE};

	// z up sensor
	auto srcArea = std::make_shared<SphereArea>(Z / 2, Real3({X / 2, Y / 2, Z}));
	real frequency = 10e+6, T = 1.0 / frequency;
	real A = -1e+6;
	borderCondition.area = srcArea;
	borderCondition.values = {
			{PhysicalQuantities::T::Szz, [A, T](real t) {return (t < T) ? A : 0; }},
			{PhysicalQuantities::T::Sxz, [] (real) {return 0; }},
			{PhysicalQuantities::T::Syz, [] (real) {return 0; }},
	};
	
	statement.cubicGridBorderConditions.push_back(borderCondition);
	
	task.statements.push_back(statement);
	return task;
}


inline Task parseTaskCgalAnisotropy() {
	Task task;

	task.modelId = Models::T::ELASTIC3D;
	task.materialId = Materials::T::ORTHOTROPIC;
	task.gridId = Grids::T::CGAL;
	task.snapshottersId = {
			Snapshotters::T::VTK
	};

	task.cgal3DGrid.spatialStep = 0.003;
//	task.cgal3DGrid.detectSharpEdges = true;
	task.cgal3DGrid.polyhedronFileName = "meshes/layers_with_fracture.off";
	
	Statement statement;
	/*real rho = 4;
	real lambda = 2;
	real mu = 1;
	statement.materialConditions.defaultMaterial =
	        std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);*/

	statement.materialConditions.defaultMaterial = createMaterial(0);
	statement.materialConditions.materials = {
			{createArea(0),  createMaterial( M_PI / 2)},
	};

	statement.globalSettings.CourantNumber = 0.5;

	statement.globalSettings.numberOfSnaps = 50;
	statement.globalSettings.stepsPerSnap = 10;

	Statement::BorderCondition borderConditionBottom;
	borderConditionBottom.area = std::make_shared<AxisAlignedBoxArea>(
			Real3({0.01, 0.01, -1}), Real3({0.15, 0.15, 0.01}));
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
	sourceUp.area = std::make_shared<StraightBoundedCylinderArea>(
			0.02, Real3({0.08, 0.08, 0.03}), Real3({0.08, 0.08, 0.05}));
	sourceUp.type = BorderConditions::T::FIXED_FORCE;
	sourceUp.values = {
		[] (real) { return 0; },
		[] (real) { return 0; },
		[] (real t) { return (t < 10e-6 /*0.005*/) ? -1 : 0; }
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
	                              sourceUp,
	                              fracture};

	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {
		PhysicalQuantities::T::PRESSURE,
		PhysicalQuantities::T::Szz,
	};

	task.statements.push_back(statement);
	return task;
}
