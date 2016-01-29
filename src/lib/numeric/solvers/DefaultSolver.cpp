#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

template<class TGrid>
void DefaultSolver<TGrid>::initialize(const Task &task) {
	LOG_INFO("Start initialization");
	this->method = new GridCharacteristicMethod<TGrid>();
	this->plasticFlowCorrector = new IdealPlasticFlowCorrector<TGrid>();
	plasticFlowCorrector->initialize(task);

	this->mesh = new TGrid();
	this->mesh->initialize(task);
	this->newMesh = new TGrid();
	this->newMesh->initialize(task);

	CourantNumber = task.CourantNumber;
	splittingSecondOrder = task.splittingSecondOrder;
}

template<class TGrid>
DefaultSolver<TGrid>::~DefaultSolver() {
	delete this->method;
	delete this->mesh;
	delete this->newMesh;
}

template<class TGrid>
void DefaultSolver<TGrid>::nextTimeStepImpl() {
	LOG_INFO("Start time step " << step);
	real tau = calculateTau();

	if (splittingSecondOrder) {
		switch (TGrid::DIMENSIONALITY) {
			case 1:
				stage(0, tau);
				break;
			case 2:
				stage(0, tau / 2);
				stage(1, tau);
				stage(0, tau / 2);
				break;
			case 3:
				THROW_UNSUPPORTED("TODO splitting second order in 3D");
				break;
		}
	} else {
		for (int s = 0; s < TGrid::DIMENSIONALITY; s++) {
			stage(s, tau);
		}
	}

	plasticFlowCorrector->apply(mesh);

}

template<class TGrid>
void DefaultSolver<TGrid>::stage(const int s, const real timeStep) {
	mesh->beforeStage();
	method->stage(s, timeStep, mesh, newMesh); // now actual values is by pointer newMesh
	std::swap(mesh, newMesh); // now actual values is again by pointer mesh
	mesh->afterStage();
}

template<class TGrid>
real DefaultSolver<TGrid>::calculateTau() const {
	return CourantNumber * mesh->getMinimalSpatialStep() / mesh->getMaximalLambda();
}

template class DefaultSolver<StructuredGrid<Elastic1DModel>>;
template class DefaultSolver<StructuredGrid<Elastic2DModel>>;
template class DefaultSolver<StructuredGrid<Elastic3DModel>>;
template class DefaultSolver<StructuredGrid<OrthotropicElastic3DModel>>;
