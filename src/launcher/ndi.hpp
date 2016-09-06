#include <libgcm/util/task/Task.hpp>

using namespace gcm;


inline real lambda(const real E, const real nu) {
	return E * nu / (1 + nu) / (1 - 2 * nu);
}

inline real mu(const real E, const real nu) {
	return E / (2 + 2 * nu);
}


inline std::shared_ptr<OrthotropicMaterial>
createCompositeMaterial(const real phi) {
	real E11 = 1e+11;
	real E22 = 3e+10;
	real nu = 0.43;
	real rho = 1.6;
	
	real c11 = lambda(E11, nu) + 2 * mu(E11, nu);
	real c22 = lambda(E22, nu) + 2 * mu(E22, nu);
	real c33 = c22;
	real c12 = lambda(E11, nu);
	real c13 = c12;
	real c23 = lambda(E22, nu);
	real c44 = mu(E22, nu);
	real c55 = mu(E11, nu);
	real c66 = c55;
	OrthotropicMaterial material(rho,
			{c11, c12, c13,
			      c22, c23,
			           c33,
			               c44, c55, c66});
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


inline Task ndi() {
	Task task;
	task.globalSettings.dimensionality = 2;
	task.globalSettings.gridId = Grids::T::CUBIC;
	task.globalSettings.snapshottersId = { Snapshotters::T::VTK };
	task.globalSettings.numberOfSnaps = 70;
	task.globalSettings.CourantNumber = 0.9;
	
	task.bodies = {
			{0, {Materials::T::ORTHOTROPIC, Models::T::ELASTIC, {}}},
			{1, {Materials::T::ISOTROPIC,   Models::T::ELASTIC, {}}}
	};
	
	task.materialConditions.type = Task::MaterialCondition::Type::BY_BODIES;
	task.materialConditions.byBodies.bodyMaterialMap = {
		{0, createCompositeMaterial(0)},
		{1, createNdiMaterial()},
	};
	
	
	task.cubicGrid.borderSize = 2;
	task.cubicGrid.h = {1, 0.25};
	task.cubicGrid.cubics = {
			{0, {{20, 41}, { 0,  0}}},
			{1, {{10, 41}, { 5, 40}}}
	};
	
	
	Task::CubicBorderCondition::Values freeBorderY = {
		{PhysicalQuantities::T::Syy, [](real) {return 0;} },
		{PhysicalQuantities::T::Sxy, [](real) {return 0;} },
	};
	Task::CubicBorderCondition::Values freeBorderX = {
		{PhysicalQuantities::T::Sxx, [](real) {return 0;} },
		{PhysicalQuantities::T::Sxy, [](real) {return 0;} },
	};
	auto area = std::make_shared<InfiniteArea>();
	task.cubicBorderConditions = {
		{0, {
			 {1, area, freeBorderY},
		 }},
		{1, {
			 {0, area, freeBorderX},
			 {1, area, freeBorderY},
		 }},
	};
	
	
	Task::InitialCondition::Wave wave;
	wave.waveType = Waves::T::P_FORWARD;
	wave.direction = 1; // along y
	wave.quantity = PhysicalQuantities::T::PRESSURE;
	wave.quantityValue = 1;
	Real3 min({-1000, 2.5, -1000});
	Real3 max({ 1000, 7.5,  1000});
	wave.area = std::make_shared<AxisAlignedBoxArea>(min, max);
	task.initialCondition.waves.push_back(wave);
	
	return task;
}

