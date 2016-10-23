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
	real m = 2.5; // times to increase acoustic impedance when acoustic velocity stay constant
	OrthotropicMaterial material(1.6e+3 * m,
			{23.3e+9 * m, 6.96e+9 * m, 6.96e+9 * m,
			         12.0e+9 * m,  6.96e+9 * m,
			                  12.0e+9 * m,
			 1.67e+9 * m, 5.01e+9 * m, 5.01e+9 * m});
	material.anglesOfRotation = {0, 0, phi};
	material.tau0 = 5e-6;
	return std::make_shared<OrthotropicMaterial>(material);
}

//inline std::shared_ptr<IsotropicMaterial>
//createLeadMetaniobate() {
///// create "lead metaniobate" material -- piezoelectric
//	real E = 68e+9; // Young
//	real nu = 0.19; // Puasson
//	real rho = 8e+3;
//	return std::make_shared<IsotropicMaterial>(rho, lambda(E, nu), mu(E, nu));
//}

inline std::shared_ptr<IsotropicMaterial>
createRexolite() {
/// create "rexolite" material -- material of NDI prism
	real lambda = 3.03e+9;
	real mu = 1.39e+9;
	real rho = 1.05e+3;
	return std::make_shared<IsotropicMaterial>(rho, lambda, mu);
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
	task.globalSettings.CourantNumber = 0.9;
	
	return task;
}


inline real ndiEmpty(Task& task, int& prismCenterX, const int nodesPerMm) {
	task = ndiCommon();
	task.globalSettings.stepsPerSnap = 5;
	task.globalSettings.numberOfSnaps = 
			75 * nodesPerMm / task.globalSettings.stepsPerSnap;
	
	task.bodies = {
			{1, {Materials::T::ISOTROPIC, Models::T::ELASTIC, {}}}
	};
	
	task.materialConditions.type = Task::MaterialCondition::Type::BY_BODIES;
	task.materialConditions.byBodies.bodyMaterialMap = {
			{1, createRexolite()},
	};
	
	real prismDiameter = 6.35e-3;
	real lengthToDiameterRatio = 1.4;
	real prismLength = prismDiameter * lengthToDiameterRatio;
	int prismSizeY = (int)(prismLength * 1e+3) * nodesPerMm;
	int prismSizeX = (int)(prismDiameter * 1e+3) * nodesPerMm;
	int prismStartX = prismSizeX / 2;
	prismCenterX = prismStartX + prismSizeX / 2;
	real hX = prismDiameter / prismSizeX;
	real hY = prismLength / prismSizeY;
	task.cubicGrid.h = {hX, hY};
	task.cubicGrid.cubics = {
			{1, {{prismSizeX, prismSizeY}, { prismStartX, 0}}}
	};
	real upperY = (prismSizeY - 1) * hY;
	real big = 1, eps = 1e-7;
//	real startX = -big, finishX = big;
	real startX = prismStartX * hX + 0.2 * prismDiameter;
	real finishX = prismStartX * hX + 0.8 * prismDiameter;
	
	
	auto upper = std::make_shared<AxisAlignedBoxArea>(
			Real3({startX, upperY - eps, -big}), Real3({finishX, upperY + eps, big}));
	real tau = 1e-7/* / 2*/;
//	real omega = 2 * M_PI / (2 * tau);
	task.cubicBorderConditions = {
		{1, {
//			 {0, std::make_shared<InfiniteArea>(), freeBorder2D(0) },
			 {1, std::make_shared<InfiniteArea>(), freeBorder2D(1) },
			 {1, upper, fixedNormalForce2D(1, [=](real t) {
					t -= 2 * tau;
					return - /*sin(omega * t - M_PI / 6) **/ exp(-t*t / ( 2 * tau*tau));}) },
		 }},
	};
	
//	Task::InitialCondition::Wave wave;
//	wave.waveType = Waves::T::P_FORWARD;
//	wave.direction = 1;
//	wave.quantity = PhysicalQuantities::T::PRESSURE;
//	wave.quantityValue = 1;
//	wave.area = std::make_shared<AxisAlignedBoxArea>(Real3({-10, 4e-3, -10}), Real3({10, 4e-3, 10}));
//	task.initialCondition.waves.push_back(wave);
	
	task.detector.area = upper;
	return hY;
}


inline Task ndi(real thickness, const int nodesPerMm) {
	Task task;
	int prismCenterX = 0;
	const real hY = ndiEmpty(task, prismCenterX, nodesPerMm);
	
	task.bodies.insert(
			{0, {Materials::T::ORTHOTROPIC, Models::T::ELASTIC, {Odes::T::MAXWELL_VISCOSITY}}});
	
	task.materialConditions.byBodies.bodyMaterialMap.insert(
			{0, createCompositeMaterial(0)});
	
	if (thickness == 0) {
		task.cubicGrid.cubics.insert({0, {{11, 11}, { 0, -111}}});
		return task;
	}
	bool nonReflection = false;
	if (thickness >= 100) {
		nonReflection = true;
		thickness = 5;
	}
	
	int sampleSizeY = (int)(thickness / 1e+3 / hY);
	assert_gt(sampleSizeY, 3);
	task.cubicGrid.cubics.insert(
			{0, {{prismCenterX * 2, sampleSizeY}, { 0, - sampleSizeY}}});
	
	std::shared_ptr<Area> area = std::make_shared<InfiniteArea>();
	if (nonReflection) {
		area = std::make_shared<AxisAlignedBoxArea>(
				Real3({-10, -hY - 1e-7, -10}), Real3({10, -hY + 1e-7, 10}));
	}
	task.cubicBorderConditions.insert({0, { {1, area, freeBorder2D(1)}, }});
	
	return task;
}


