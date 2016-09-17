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
	real E = 68e+9/2; // Young
	real nu = 0.19; // Puasson
	real rho = 8e+3;
	return std::make_shared<IsotropicMaterial>(rho, lambda(E, nu), mu(E, nu));
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
	task.globalSettings.stepsPerSnap = 7;
	
	task.bodies = {
			{1, {Materials::T::ISOTROPIC,   Models::T::ELASTIC, {}}}
	};
	
	task.materialConditions.type = Task::MaterialCondition::Type::BY_BODIES;
	task.materialConditions.byBodies.bodyMaterialMap = {
		{1, createNdiMaterial()},
	};
	
	real prismDiameter = 6.35e-3;
	real prismLength = 6e-3;
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
					IsotropicMaterial(1e+3, 12e+9, 3e+9)) });
	
	int sampleSizeY = 41;
	task.cubicGrid.cubics.insert(
			{0, {{21, sampleSizeY}, { 0, 1 - sampleSizeY - 100}}});
	
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


