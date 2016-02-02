#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

template<class TGrid>
void DefaultSolver<TGrid>::initializeImpl(const Task &task) {
	LOG_INFO("Start initialization");
	method = new GridCharacteristicMethod<TGrid>();
	corrector = new typename TGrid::Model::Corrector();
	corrector->initialize(task);
	internalOde = new typename TGrid::Model::InternalOde();
	internalOde->initialize(task);

	mesh = new TGrid();
	mesh->initialize(task);
	newMesh = new TGrid();
	newMesh->initialize(task);

	CourantNumber = task.CourantNumber;
	splittingSecondOrder = task.splittingSecondOrder;
}

template<class TGrid>
DefaultSolver<TGrid>::~DefaultSolver() {
	delete method;
	delete corrector;
	delete internalOde;
	delete mesh;
	delete newMesh;
}

template<class TGrid>
void DefaultSolver<TGrid>::nextTimeStepImpl() {
	LOG_INFO("Start time step " << step);
	mesh->beforeStep();
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
	fixVariablesOrder();

	internalOdeNextStep(tau);
	applyCorrectors();
	mesh->afterStep();
}

template<class TGrid>
void DefaultSolver<TGrid>::stage(const int s, const real timeStep) {
	mesh->beforeStage();
	method->stage(s, timeStep, mesh, newMesh); // now actual PDE values is by pointer newMesh
	odeShiftedFromPde = !odeShiftedFromPde; // but actual ODE values wasn't moved
	std::swap(mesh, newMesh); // now actual PDE values is again by pointer mesh
	mesh->afterStage();
}

template<class TGrid>
void DefaultSolver<TGrid>::fixVariablesOrder() {
	if (odeShiftedFromPde /*&& TGrid::Model::InternalOde::NonTrivial*/) {
		assert_eq(mesh->nodes.size(), newMesh->nodes.size());
		for (unsigned long i = 0; i < mesh->nodes.size(); i++) {
			std::swap(mesh->nodes[i].u, newMesh->nodes[i].u);
		}
		std::swap(mesh, newMesh);
	}
}

template<class TGrid>
void DefaultSolver<TGrid>::internalOdeNextStep(const real timeStep) {
	if (TGrid::Model::InternalOde::NonTrivial) {
		assert_eq(mesh->nodes.size(), newMesh->nodes.size());
		for (unsigned long i = 0; i < mesh->nodes.size(); i++) {
			internalOde->nextStep(mesh->nodes[i], newMesh->nodes[i], timeStep);
		}
	}
};

template<class TGrid>
void DefaultSolver<TGrid>::applyCorrectors() {
	if (TGrid::Model::Corrector::NonTrivial) {
		for (auto& node : mesh->nodes) {
			corrector->apply(node);
		}
	}
}

template<class TGrid>
real DefaultSolver<TGrid>::calculateTau() const {
	return CourantNumber * mesh->getMinimalSpatialStep() / mesh->getMaximalLambda();
}

template class DefaultSolver<StructuredGrid<Elastic1DModel>>;
template class DefaultSolver<StructuredGrid<Elastic2DModel>>;
template class DefaultSolver<StructuredGrid<Elastic3DModel>>;
template class DefaultSolver<StructuredGrid<OrthotropicElastic3DModel>>;
template class DefaultSolver<StructuredGrid<ContinualDamageElastic2DModel>>;
template class DefaultSolver<StructuredGrid<IdealPlastic2DModel>>;

template class DefaultSolver<StructuredGrid<SuperDuperModel>>;

