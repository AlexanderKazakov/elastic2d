#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

template<class TMesh>
void DefaultSolver<TMesh>::initializeImpl(const Task &task) {
	LOG_INFO("Start initialization");

	mesh = new TMesh();
	mesh->initialize(task);

	borderConditions = new Border();
	borderConditions->initialize(task);

	splittingSecondOrder = task.splittingSecondOrder;
}

template<class TMesh>
DefaultSolver<TMesh>::~DefaultSolver() {
	delete borderConditions;
	delete mesh;
}

template<class TMesh>
void DefaultSolver<TMesh>::beforeStatementImpl(const Statement& statement) {
	CourantNumber = statement.CourantNumber;
	corrector = new Corrector();
	corrector->beforeStatement(statement);
	internalOde = new InternalOde();
	internalOde->beforeStatement(statement);

	borderConditions->beforeStatement(statement);
	mesh->beforeStatement(statement);
}

template<class TMesh>
void DefaultSolver<TMesh>::nextTimeStepImpl() {
	LOG_INFO("Start time step " << step);
	mesh->beforeStep();
	real tau = calculateTau();

	if (splittingSecondOrder) {
		switch (TMesh::DIMENSIONALITY) {
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
		for (int s = 0; s < TMesh::DIMENSIONALITY; s++) {
			stage(s, tau);
		}
	}

	internalOdeNextStep(tau);
	applyCorrectors();
	moveMesh(tau);
	mesh->afterStep();
}

template<class TMesh>
void DefaultSolver<TMesh>::afterStatementImpl() {
	mesh->afterStatement();
	delete corrector;
	delete internalOde;
}

template<class TMesh>
void DefaultSolver<TMesh>::stage(const int s, const real timeStep) {
	DataBus<Model, Grid>::exchangeNodesWithNeighbors(mesh);
	borderConditions->applyBorderBeforeStage(mesh, getCurrentTime(), timeStep, s);
	GCM::stage(s, timeStep, mesh); // now actual PDE values is in pdeVectorsNew
	borderConditions->applyBorderAfterStage(mesh, getCurrentTime(), timeStep, s);
	std::swap(mesh->pdeVectors, mesh->pdeVectorsNew); // return actual PDE values back to pdeVectors
}

template<class TMesh>
void DefaultSolver<TMesh>::internalOdeNextStep(const real timeStep) {
	if (InternalOde::NonTrivial) {
		assert_eq(mesh->pdeVectors.size(), mesh->odeValues.size());
		for (auto it : *mesh) {
			internalOde->nextStep(mesh->_ode(it), mesh->pde(it), timeStep);
		}
	}
}

template<class TMesh>
void DefaultSolver<TMesh>::applyCorrectors() {
	if (Corrector::NonTrivial) {
		for (auto it : *mesh) {
			corrector->apply(mesh->_pde(it));
		}
	}
}

template<class TMesh>
void DefaultSolver<TMesh>::moveMesh(const real timeStep) {
	MeshMover<Model, Grid>::moveMesh(*mesh, timeStep);
}

template<class TMesh>
real DefaultSolver<TMesh>::calculateTau() const {
	return CourantNumber * mesh->getMinimalSpatialStep() / mesh->getMaximalLambda();
}

template class DefaultSolver<DefaultMesh<Elastic1DModel, CubicGrid>>;
template class DefaultSolver<DefaultMesh<Elastic2DModel, CubicGrid>>;
template class DefaultSolver<DefaultMesh<Elastic3DModel, CubicGrid>>;
template class DefaultSolver<DefaultMesh<OrthotropicElastic3DModel, CubicGrid>>;
template class DefaultSolver<DefaultMesh<ContinualDamageElastic2DModel, CubicGrid>>;
template class DefaultSolver<DefaultMesh<IdealPlastic2DModel, CubicGrid>>;
template class DefaultSolver<DefaultMesh<SuperDuperModel, CubicGrid>>;

template class DefaultSolver<DefaultMesh<Elastic2DModel, Cgal2DGrid>>;

