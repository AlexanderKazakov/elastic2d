#ifndef LIBGCM_DEFAULTSOLVER_HPP
#define LIBGCM_DEFAULTSOLVER_HPP

#include <lib/numeric/solvers/Solver.hpp>
#include <lib/numeric/gcm/GridCharacteristicMethod.hpp>
#include <lib/numeric/border_conditions/BorderConditions.hpp>
#include <lib/util/Logging.hpp>
#include <lib/util/task/Task.hpp>
#include <lib/mesh/grid/AbstractGrid.hpp>
#include <lib/mesh/DataBus.hpp>
#include <lib/mesh/MeshMover.hpp>
#include <lib/rheology/correctors/correctors.hpp>

namespace gcm {
/**
 * Class for handling complete time step
 * @tparam TMesh type of mesh to deal with
 */
template<class TMesh>
class DefaultSolver : public Solver {
public:
	typedef TMesh                                   Mesh;
	typedef typename Mesh::Model                    Model;
	typedef typename Mesh::Grid                     Grid;
	typedef typename Mesh::Material                 Material;

	typedef typename Model::Corrector               Corrector;
	typedef typename Model::InternalOde             InternalOde;

	typedef OldBorderConditions<Model, Grid, Material> Border;
	typedef DataBus<Model, Grid, Material>          DATA_BUS;
	typedef MeshMover<Model, Grid, Material>        MESH_MOVER;
	typedef GridCharacteristicMethod<Mesh>          GCM;

	DefaultSolver(const Task& task);
	virtual ~DefaultSolver();

	virtual void beforeStatement(const Statement& statement) override;
	virtual void afterStatement() override;

	virtual void nextTimeStep() override;

	/** @return mesh with actual values */
	virtual AbstractGrid* getActualMesh() const { return mesh; }

	/** Calculate time step from Courant–Friedrichs–Lewy condition */
	virtual real calculateTimeStep() const override;

protected:
	real CourantNumber = 0.0; ///< number from Courant–Friedrichs–Lewy condition

	Corrector* corrector = nullptr;
	InternalOde* internalOde = nullptr;
	Border* borderConditions = nullptr;

	Mesh* mesh = nullptr;

	void stage(const int s, const real timeStep);
	void internalOdeNextStep(const real timeStep);
	void applyCorrectors();
	void moveMesh(const real timeStep);

	USE_AND_INIT_LOGGER("gcm.DefaultSolver")
};


template<class TMesh>
DefaultSolver<TMesh>::
DefaultSolver(const Task& task) : Solver(task) {
	mesh = new Mesh(task);
	borderConditions = new Border(task);
}


template<class TMesh>
DefaultSolver<TMesh>::~DefaultSolver() {
	delete borderConditions;
	delete mesh;
}


template<class TMesh>
void DefaultSolver<TMesh>::
beforeStatement(const Statement& statement) {
	CourantNumber = statement.globalSettings.CourantNumber;

	corrector = new Corrector(statement);
	internalOde = new InternalOde(statement);

	borderConditions->beforeStatement(statement);
	mesh->beforeStatement(statement);
}


template<class TMesh>
void DefaultSolver<TMesh>::
nextTimeStep() {
	real timeStep = Clock::TimeStep();
	for (int s = 0; s < Mesh::DIMENSIONALITY; s++) {
		stage(s, timeStep);
	}
	internalOdeNextStep(timeStep);
	moveMesh(timeStep);
	applyCorrectors();
}


template<class TMesh>
void DefaultSolver<TMesh>::
afterStatement() {
	mesh->afterStatement();
	delete corrector;
	delete internalOde;
}


template<class TMesh>
void DefaultSolver<TMesh>::
stage(const int s, const real timeStep) {
	LOG_DEBUG("Start stage " << s << " ... ");
	DATA_BUS::exchangeNodesWithNeighbors(mesh);
	
	borderConditions->applyBorderBeforeStage(mesh, timeStep, s);
	GCM::stage(s, timeStep, mesh); // now actual PDE values is in pdeVectorsNew
	borderConditions->applyBorderAfterStage(mesh, timeStep, s);
	
	std::swap(mesh->pdeVectors, mesh->pdeVectorsNew); // return actual PDE values back to
		                                          // pdeVectors
}


template<class TMesh>
void DefaultSolver<TMesh>::
internalOdeNextStep(const real timeStep) {
	if (InternalOde::NonTrivial) {
		assert_eq(mesh->pdeVectors.size(), mesh->odeValues.size());
		for (auto it : *mesh) {
			internalOde->nextStep(mesh->node(it), timeStep);
		}
	}
}


template<class TMesh>
void DefaultSolver<TMesh>::
applyCorrectors() {
	if (Corrector::NonTrivial) {
		for (auto it : *mesh) {
			corrector->apply(mesh->node(it));
		}
	}
}


template<class TMesh>
void DefaultSolver<TMesh>::
moveMesh(const real timeStep) {
	MESH_MOVER::moveMesh(*mesh, timeStep);
}


template<class TMesh>
real DefaultSolver<TMesh>::
calculateTimeStep() const {
	return CourantNumber * mesh->getMinimalSpatialStep() / mesh->getMaximalEigenvalue();
}


}

#endif // LIBGCM_DEFAULTSOLVER_HPP
