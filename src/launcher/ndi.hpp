#include <libgcm/util/task/Task.hpp>

using namespace gcm;


inline Task::CubicBorderCondition::Values
fixedNormalForce2D(
		const int direction, const Task::TimeDependency normalForce) {
	switch (direction) {
		case 0:
			return {
				{PhysicalQuantities::T::Sxx, normalForce },
				{PhysicalQuantities::T::Sxy, [](real) {return 0;} },
			};
		case 1:
			return {
				{PhysicalQuantities::T::Syy, normalForce },
				{PhysicalQuantities::T::Sxy, [](real) {return 0;} },
			};
		default:
			THROW_INVALID_ARG("Invalid direction value");
	}
}

inline Task::CubicBorderCondition::Values
freeBorder2D(const int direction) {
	return fixedNormalForce2D(direction, [](real) {return 0;});
}

inline Task::CubicBorderCondition::Values
fixedNormalForce3D(
		const int direction, const Task::TimeDependency normalForce) {
	switch (direction) {
		case 0:
			return {
				{PhysicalQuantities::T::Sxx, normalForce },
				{PhysicalQuantities::T::Sxy, [](real) {return 0;} },
				{PhysicalQuantities::T::Sxz, [](real) {return 0;} },
			};
		case 1:
			return {
				{PhysicalQuantities::T::Syy, normalForce },
				{PhysicalQuantities::T::Sxy, [](real) {return 0;} },
				{PhysicalQuantities::T::Syz, [](real) {return 0;} },
			};
		case 2:
			return {
				{PhysicalQuantities::T::Szz, normalForce },
				{PhysicalQuantities::T::Sxz, [](real) {return 0;} },
				{PhysicalQuantities::T::Syz, [](real) {return 0;} },
			};
		default:
			THROW_INVALID_ARG("Invalid direction value");
	}
}

inline Task::CubicBorderCondition::Values
freeBorder3D(const int direction) {
	return fixedNormalForce3D(direction, [](real) {return 0;});
}

inline Task::CubicBorderCondition::Values
fixedNormalVelocity3D(
		const int direction, const Task::TimeDependency normalVelocity) {
	switch (direction) {
		case 0:
			return {
				{PhysicalQuantities::T::Vx, normalVelocity},
				{PhysicalQuantities::T::Vy, [](real) {return 0;} },
				{PhysicalQuantities::T::Vz, [](real) {return 0;} },
			};
		case 1:
			return {
				{PhysicalQuantities::T::Vx, [](real) {return 0;} },
				{PhysicalQuantities::T::Vy, normalVelocity},
				{PhysicalQuantities::T::Vz, [](real) {return 0;} },
			};
		case 2:
			return {
				{PhysicalQuantities::T::Vx, [](real) {return 0;} },
				{PhysicalQuantities::T::Vy, [](real) {return 0;} },
				{PhysicalQuantities::T::Vz, normalVelocity},
			};
		default:
			THROW_INVALID_ARG("Invalid direction value");
	}
}


inline real lambda(const real E, const real nu) {
	return E * nu / (1 + nu) / (1 - 2 * nu);
}

inline real mu(const real E, const real nu) {
	return E / (2 + 2 * nu);
}


inline std::shared_ptr<OrthotropicMaterial>
createCompositeMaterial(const real phi) {
	real rho = 1.6e+3;
	OrthotropicMaterial material(rho,
			{23.3e+9, 6.96e+9, 6.96e+9,
			         10.3e+9,  6.96e+9,
			                  10.3e+9, 1.67, 5.01, 5.01});
	material.anglesOfRotation = {0, 0, phi};
	return std::make_shared<OrthotropicMaterial>(material);
}

inline std::shared_ptr<IsotropicMaterial>
createNdiMaterial() {
/// create "lead metaniobate" material for ndi prism
	real E = 68e+9; // Young
	real nu = 0.19; // Puasson
	real rho = 6e+3;
	return std::make_shared<IsotropicMaterial>(rho, lambda(E, nu), mu(E, nu));
}


inline std::shared_ptr<OrthotropicMaterial>
createCompositeMaterialFromCagiReport2014() {
	real rho = 1.58e+3;
	OrthotropicMaterial material(rho,
			{10.30e+9,  6.96e+9,  6.96e+9,
			           23.25e+9,  6.96e+9,
			                     10.30e+9,
			                               5.01e+9, 1.67e+9, 5.01e+9});
	material.anglesOfRotation = {0, 0, 0};
	material.materialNumber = 0;
	return std::make_shared<OrthotropicMaterial>(material);
}

inline std::shared_ptr<IsotropicMaterial>
createTitan() {
	real E = 120e+9; // Young
	real nu = 0.31; // Puasson
	real rho = 4.5e+3;
	std::shared_ptr<IsotropicMaterial> ans =
			std::make_shared<IsotropicMaterial>(rho, lambda(E, nu), mu(E, nu));
	ans->materialNumber = 1;
	return ans;
}

inline std::shared_ptr<IsotropicMaterial>
createSteel() {
	real lambda = 99.43e+9;
	real mu = 78.13e+9;
	real rho = 7.8e+3;
	std::shared_ptr<IsotropicMaterial> ans =
			std::make_shared<IsotropicMaterial>(rho, lambda, mu);
	ans->materialNumber = 2;
	return ans;
}

inline std::shared_ptr<OrthotropicMaterial>
createTitanOrthotropicMaterial() {
/// isotropic is orthotropic too
	return std::make_shared<OrthotropicMaterial>(*createTitan());
}


inline Task ndiCommon() {
	Task task;
	task.cubicGrid.borderSize = 2;
	task.globalSettings.dimensionality = 2;
	task.globalSettings.gridId = Grids::T::CUBIC;
	task.globalSettings.snapshottersId = {
			Snapshotters::T::VTK,
			Snapshotters::T::SLICESNAP
	};
	task.vtkSnapshotter.quantitiesToSnap = { PhysicalQuantities::T::PRESSURE };
	task.detector.quantities = { PhysicalQuantities::T::Vy };
	task.detector.gridId = 1;
	task.globalSettings.CourantNumber = 0.5;
	
	return task;
}


inline Task ndiEmpty() {
	Task task = ndiCommon();
	task.globalSettings.numberOfSnaps = 100;
	task.globalSettings.stepsPerSnap = 5;
	
	task.bodies = {
			{1, {Materials::T::ISOTROPIC,   Models::T::ELASTIC, {}}}
	};
	
	task.materialConditions.type = Task::MaterialCondition::Type::BY_BODIES;
	task.materialConditions.byBodies.bodyMaterialMap = {
		{1, createNdiMaterial()},
	};
	
	real prismDiameter = 6.35e-3;
	real prismLength = 10e-3;
	int prismSizeY = 41, prismSizeX = 11;
	real hX = prismDiameter / prismSizeX, hY = prismLength / prismSizeY;
	task.cubicGrid.h = {hX, hY};
	task.cubicGrid.cubics = {
			{1, {{prismSizeX, prismSizeY}, { 5, 0}}}
	};
	real upperY = (prismSizeY - 1) * hY;
	real big = 1e+5, eps = 1e-5;
	
	
	auto upper = std::make_shared<AxisAlignedBoxArea>(
			Real3({-big, upperY - eps, -big}), Real3({big, upperY + eps, big}));
	real tau = 0.55e-6 / 4;
//	real omega = 2 * M_PI * 3.66e+6;
	task.cubicBorderConditions = {
		{1, {
			 {0, std::make_shared<InfiniteArea>(), freeBorder2D(0) },
			 {1, std::make_shared<InfiniteArea>(), freeBorder2D(1) },
			 {1, upper, fixedNormalForce2D(1, [=](real t) {
					t -= 2 * tau;
					return /*sin(omega * t) **/ -exp(-t*t / ( 2 * tau*tau));}) },
		 }},
	};
	
	task.detector.area = upper;
	return task;
}


inline Task ndi() {
	Task task = ndiEmpty();
	
	task.bodies.insert(
			{0, {Materials::T::ORTHOTROPIC, Models::T::ELASTIC, {}}});
	
	task.materialConditions.byBodies.bodyMaterialMap.insert(
//			{0, createCompositeMaterial(0)});
			{0, std::make_shared<OrthotropicMaterial>(
					IsotropicMaterial(1.6e+3, 7e+9, 1.5e+9)) });
	
	int sampleSizeY = 31;
	task.cubicGrid.cubics.insert(
			{0, {{21, sampleSizeY}, { 0, 1 - sampleSizeY}}});
	
	task.cubicBorderConditions.insert(
		{0, {
			 {1, std::make_shared<InfiniteArea>(), freeBorder2D(1)},
		 }});
	
	
//	Task::InitialCondition::Quantity pressure;
//	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
//	pressure.value = 0.5;
//	pressure.area = std::make_shared<SphereArea>(1e-3, Real3({6e-3, -3e-3, 0}));
//	task.initialCondition.quantities.push_back(pressure);
	
	return task;
}


inline Task titan() {
	Task task;
	task.cubicGrid.borderSize = 2;
	task.globalSettings.dimensionality = 3;
	task.globalSettings.gridId = Grids::T::CUBIC;
	task.bodies = {
			{1, {Materials::T::ORTHOTROPIC, Models::T::ELASTIC, {}}}
	};
	
	task.globalSettings.snapshottersId = {
			Snapshotters::T::VTK,
//			Snapshotters::T::SLICESNAP
	};
	task.vtkSnapshotter.quantitiesToSnap = {
			PhysicalQuantities::T::PRESSURE,
			PhysicalQuantities::T::Sxx,
			PhysicalQuantities::T::Sxy,
			PhysicalQuantities::T::Sxz,
			PhysicalQuantities::T::Syy,
			PhysicalQuantities::T::Syz,
			PhysicalQuantities::T::Szz
	};
//	task.detector.quantities = { PhysicalQuantities::T::Vy };
//	task.detector.gridId = 1;
	
	task.globalSettings.CourantNumber = 1;
	task.globalSettings.numberOfSnaps = 50;
	task.globalSettings.stepsPerSnap = 5;
	
	const real width = 4e-2;
	const real compositeH = 6.5e-3;
	const real titanH = 1e-3;
	task.materialConditions.type = Task::MaterialCondition::Type::BY_AREAS;
	task.materialConditions.byAreas.defaultMaterial =
			createCompositeMaterialFromCagiReport2014();
	Task::MaterialCondition::ByAreas::Inhomogenity titanArea = {
		std::make_shared<AxisAlignedBoxArea>(Real3({-width, 0, -width}), Real3({width, titanH + 1e-3, width})),
		createTitanOrthotropicMaterial()
	};
	task.materialConditions.byAreas.materials = {titanArea};
	
	const int sizeXZ = 50;
	const int sizeYofComposite = 50;
	const int sizeY = int(sizeYofComposite * (compositeH + titanH) / compositeH);
	const real hXZ = width / sizeXZ;
	const real hY = compositeH / sizeYofComposite;
	task.cubicGrid.h = {hXZ, hY, hXZ};
	task.cubicGrid.cubics = {
			{1, {{sizeXZ, sizeY, sizeXZ}, { -sizeXZ/2, -sizeYofComposite, -sizeXZ/2}}}
	};
	
	const real diameter = 3e-3;
	auto upper = std::make_shared<StraightBoundedCylinderArea>(diameter / 2,
			Real3({0, -compositeH/2, 0}), Real3({0, titanH + 1e-3, 0}));
	task.cubicBorderConditions = {
		{1, {
			 {1, std::make_shared<InfiniteArea>(), freeBorder3D(1) },
			 {1, upper, fixedNormalVelocity3D(1, [=](real) {return -4.24;}) },
		 }},
	};
//	task.detector.area = upper;
	
	return task;
}

