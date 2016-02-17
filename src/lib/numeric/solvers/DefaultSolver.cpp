#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/grid/DefaultGrid.hpp>
#include <lib/grid/Cgal2DGrid.hpp>
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

	grid = new TGrid();
	grid->initialize(task);

	CourantNumber = task.CourantNumber;
	splittingSecondOrder = task.splittingSecondOrder;
}

template<class TGrid>
DefaultSolver<TGrid>::~DefaultSolver() {
	delete method;
	delete corrector;
	delete internalOde;
	delete grid;
}

template<class TGrid>
void DefaultSolver<TGrid>::nextTimeStepImpl() {
	LOG_INFO("Start time step " << step);
	grid->beforeStep();
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

	internalOdeNextStep(tau);
	applyCorrectors();
	grid->afterStep();
}

template<class TGrid>
void DefaultSolver<TGrid>::stage(const int s, const real timeStep) {
	grid->beforeStage();
	method->stage(s, timeStep, grid); // now actual PDE values is in pdeVectorsNew
	std::swap(grid->pdeVectors, grid->pdeVectorsNew); // return them back to pdeVectors
	grid->afterStage();
}

template<class TGrid>
void DefaultSolver<TGrid>::internalOdeNextStep(const real timeStep) {
	if (TGrid::Model::InternalOde::NonTrivial) {
		assert_eq(grid->pdeVectors.size(), grid->odeValues.size());
		for (auto it : *grid) {
			internalOde->nextStep(grid->_ode(it), grid->pde(it), timeStep);
		}
	}
}

template<class TGrid>
void DefaultSolver<TGrid>::applyCorrectors() {
	if (TGrid::Model::Corrector::NonTrivial) {
		for (auto it : *grid) {
			corrector->apply(grid->_pde(it));
		}
	}
}

template<class TGrid>
real DefaultSolver<TGrid>::calculateTau() const {
	return CourantNumber * grid->getMinimalSpatialStep() / grid->getMaximalLambda();
}

template class DefaultSolver<DefaultGrid<Elastic1DModel, StructuredGrid>>;
template class DefaultSolver<DefaultGrid<Elastic2DModel, StructuredGrid>>;
template class DefaultSolver<DefaultGrid<Elastic3DModel, StructuredGrid>>;
template class DefaultSolver<DefaultGrid<OrthotropicElastic3DModel, StructuredGrid>>;
template class DefaultSolver<DefaultGrid<ContinualDamageElastic2DModel, StructuredGrid>>;
template class DefaultSolver<DefaultGrid<IdealPlastic2DModel, StructuredGrid>>;

template class DefaultSolver<DefaultGrid<SuperDuperModel, StructuredGrid>>;

//template class DefaultSolver<Cgal2DGrid<Elastic2DModel>>;

