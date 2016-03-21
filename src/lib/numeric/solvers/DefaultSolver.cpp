#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

template<class TMesh>
DefaultSolver<TMesh>::DefaultSolver(const Task &task) : Solver(task) {
	LOG_INFO("Start initialization");

	mesh = new Mesh(task);
	borderConditions = new Border(task);
}

template<class TMesh>
DefaultSolver<TMesh>::~DefaultSolver() {
	delete borderConditions;
	delete mesh;
}

template<class TMesh>
void DefaultSolver<TMesh>::beforeStatement(const Statement& statement) {
	CourantNumber = statement.globalSettings.CourantNumber;
	splittingSecondOrder = statement.globalSettings.splittingSecondOrder;

	corrector = new Corrector(statement);
	internalOde = new InternalOde(statement);

	borderConditions->beforeStatement(statement);
	mesh->beforeStatement(statement);
}

template<class TMesh>
void DefaultSolver<TMesh>::nextTimeStep(const real timeStep) {
	LOG_INFO("Start time step " << step);

	for (int s = 0; s < Mesh::DIMENSIONALITY; s++) {
		stage(s, timeStep);
	}
	internalOdeNextStep(timeStep);
	moveMesh(timeStep);
	applyCorrectors();
}

template<class TMesh>
void DefaultSolver<TMesh>::afterStatement() {
	mesh->afterStatement();
	delete corrector;
	delete internalOde;
}

template<class TMesh>
void DefaultSolver<TMesh>::stage(const int s, const real timeStep) {
	DATA_BUS::exchangeNodesWithNeighbors(mesh);
	borderConditions->applyBorderBeforeStage(mesh, timeStep, s);
	GCM::stage(s, timeStep, mesh); // now actual PDE values is in pdeVectorsNew
	borderConditions->applyBorderAfterStage(mesh, timeStep, s);
	std::swap(mesh->pdeVectors, mesh->pdeVectorsNew); // return actual PDE values back to pdeVectors
}

template<class TMesh>
void DefaultSolver<TMesh>::internalOdeNextStep(const real timeStep) {
	if (InternalOde::NonTrivial) {
		assert_eq(mesh->pdeVectors.size(), mesh->odeValues.size());
		for (auto it : *mesh) {
			internalOde->nextStep(mesh->node(it), timeStep);
		}
	}
}

template<class TMesh>
void DefaultSolver<TMesh>::applyCorrectors() {
	if (Corrector::NonTrivial) {
		for (auto it : *mesh) {
			corrector->apply(mesh->node(it));
		}
	}
}

template<class TMesh>
void DefaultSolver<TMesh>::moveMesh(const real timeStep) {
	MESH_MOVER::moveMesh(*mesh, timeStep);
}

template<class TMesh>
real DefaultSolver<TMesh>::calculateTimeStep() const {
	return CourantNumber * mesh->getMinimalSpatialStep() / mesh->getMaximalEigenvalue();
}

template class DefaultSolver<DefaultMesh<Elastic1DModel, CubicGrid, IsotropicMaterial>>;
template class DefaultSolver<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>>;
template class DefaultSolver<DefaultMesh<Elastic3DModel, CubicGrid, IsotropicMaterial>>;

template class DefaultSolver<DefaultMesh<SuperDuperModel, CubicGrid, IsotropicMaterial>>;
template class DefaultSolver<DefaultMesh<SuperDuperModel, CubicGrid, OrthotropicMaterial>>;

template class DefaultSolver<DefaultMesh<Elastic2DModel, Cgal2DGrid, IsotropicMaterial>>;

