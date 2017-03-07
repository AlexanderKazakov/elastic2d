#include <libgcm/util/task/Task.hpp>
#include <libgcm/util/math/Area.hpp>

using namespace gcm;

inline Task torAcoustic() {
    Task task;
    
    task.globalSettings.dimensionality = 3;
    task.globalSettings.gridId = Grids::T::SIMPLEX;
    task.globalSettings.snapshottersId = { Snapshotters::T::VTK };
    task.globalSettings.CourantNumber = 1;
    task.globalSettings.numberOfSnaps = 100;
    task.globalSettings.stepsPerSnap = 1;
    
    task.simplexGrid.mesher = Task::SimplexGrid::Mesher::INM_MESHER;
    task.simplexGrid.fileName = "meshes/coarse/tor.out";
    task.simplexGrid.scale = 1;
    
    task.materialConditions.type = Task::MaterialCondition::Type::BY_BODIES;
    task.bodies = {
        {1, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {}}},
    };
    auto tor = std::make_shared<IsotropicMaterial>(1.904,  7.854, 0, 0, 0, 4, 0);//form bones
    task.materialConditions.byBodies.bodyMaterialMap = {
        {1, tor},
    };
    
    //Task::InitialCondition::Quantity pressure;
    //pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
    //pressure.value = 1;
    //pressure.area = std::make_shared<SphereArea>(2, Real3({0, 5, 147}));
    //task.initialCondition.quantities.push_back(pressure);
    
    task.vtkSnapshotter.quantitiesToSnap = { PhysicalQuantities::T::PRESSURE };
    
    //task.vtkSnapshotter.quantitiesToSnap = { //May I use it?
    //    PhysicalQuantities::T::PRESSURE,
    //    PhysicalQuantities::T::Sxx,
    //    PhysicalQuantities::T::Sxy,
    //    PhysicalQuantities::T::Sxz,
    //    PhysicalQuantities::T::Syy,
    //    PhysicalQuantities::T::Syz,
    //    PhysicalQuantities::T::Szz
    //};
	
	task.contactCondition.defaultCondition = ContactConditions::T::SLIDE;
	
	Task::BorderCondition freeBorder;
	freeBorder.area = std::make_shared<InfiniteArea>();
	freeBorder.type = BorderConditions::T::FIXED_FORCE;
    freeBorder.values = {
        [] (real) { return 0; },
    };
    
	Task::BorderCondition source;
	source.area = std::make_shared<SphereArea>(1, Real3({-6.5, 0, 0}));
	source.type = BorderConditions::T::FIXED_FORCE;
    source.values = {
        [] (real t) { return (t < 0.2) ? 1 : 0; },
    };
	task.borderConditions = {freeBorder, source};
	
	return task;
}