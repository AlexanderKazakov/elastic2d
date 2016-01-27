#include <lib/numeric/solvers/DefaultSolver.hpp>

using namespace gcm;

template<class TNode>
void DefaultSolver<TNode>::initialize(const Task &task) {
	LOG_INFO("Start initialization");
	this->method = new GridCharacteristicMethod<TNode>();

	this->mesh = new StructuredGrid<TNode>();
	this->mesh->initialize(task);
	this->newMesh = new StructuredGrid<TNode>();
	this->newMesh->initialize(task);

	CourantNumber = task.CourantNumber;
	splittingSecondOrder = task.splittingSecondOrder;
}

template<class TNode>
DefaultSolver<TNode>::~DefaultSolver() {
	delete this->method;
	delete this->mesh;
	delete this->newMesh;
}

template<class TNode>
void DefaultSolver<TNode>::nextTimeStep() {
	LOG_INFO("Start time step " << step);
	real tau = calculateTau();

	if (splittingSecondOrder) {
		switch (TNode::DIMENSIONALITY) {
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
		for (int s = 0; s < TNode::DIMENSIONALITY; s++) {
			stage(s, tau);
		}
	}

	currentTime += tau;
}

template<class TNode>
void DefaultSolver<TNode>::stage(const int s, const real &timeStep) {
	mesh->beforeStage();
	method->stage(s, timeStep, mesh, newMesh); // now actual values is by pointer newMesh
	std::swap(mesh, newMesh); // now actual values is again by pointer mesh
	mesh->afterStage();
}

template<class TNode>
real DefaultSolver<TNode>::calculateTau() const {
	return CourantNumber * mesh->getMinimalSpatialStep() / mesh->getMaximalLambda();
}

template class DefaultSolver<IdealElastic1DNode>;
template class DefaultSolver<IdealElastic2DNode>;
template class DefaultSolver<IdealElastic3DNode>;