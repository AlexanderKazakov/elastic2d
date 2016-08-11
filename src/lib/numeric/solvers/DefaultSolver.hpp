#ifndef LIBGCM_DEFAULTSOLVER_HPP
#define LIBGCM_DEFAULTSOLVER_HPP

#include <lib/numeric/solvers/Solver.hpp>
#include <lib/numeric/gcm/grid_characteristic_methods.hpp>
#include <lib/numeric/special_border_conditions/SpecialBorderConditions.hpp>
#include <lib/util/Logging.hpp>
#include <lib/util/task/Task.hpp>
#include <lib/mesh/grid/AbstractGrid.hpp>
#include <lib/mesh/DataBus.hpp>
#include <lib/mesh/MeshMover.hpp>
#include <lib/rheology/correctors/correctors.hpp>


namespace gcm {

/**
 * Class for handling mesh of a concrete type
 */
template<class TMesh>
class DefaultSolver : public Solver {
public:
	typedef TMesh                                           Mesh;
	typedef typename Mesh::Model                            Model;
	typedef typename Mesh::Grid                             Grid;
	typedef typename Mesh::Material                         Material;
	typedef typename Mesh::GridId                           GridId;
	typedef typename Mesh::GlobalScene                      GlobalScene;
	
	typedef typename Model::Corrector                       Corrector;
	typedef typename Model::InternalOde                     InternalOde;
	
	typedef SpecialBorderConditions<Model, Grid, Material>  SpecialBorder;
	typedef DataBus<Model, Grid, Material>                  DATA_BUS;
	typedef MeshMover<Model, Grid, Material>                MESH_MOVER;
	typedef GridCharacteristicMethod<Model, Grid, Material> GcmMethod;
	
	static const int GRID_DIMENSIONALITY = Mesh::GRID_DIMENSIONALITY;
	
	DefaultSolver(const Task& task,
			AbstractGlobalScene* abstractGlobalScene, const GridId gridId_);
	virtual ~DefaultSolver();
	
	virtual void beforeStatement(const Statement& statement) override;
	virtual void afterStatement() override;
	
	
	/**
	 * All necessary solver actions before stages would performed
	 */
	virtual void beforeStage() override;
	
	/**
	 * Calculate contact nodes in given stage
	 */
	virtual void contactStage(const int s, const real timeStep) override;
	
	/**
	 * Calculate inner and border nodes in given stage
	 */
	virtual void privateStage(const int s, const real timeStep) override;
	
	/**
	 * All necessary solver actions after stages performed
	 */
	virtual void afterStages(const real timeStep) override;
	
	
	/// abstract interface to mesh
	virtual AbstractGrid* getAbstractMesh() const override { return mesh; }
	
	/** Calculate time step from Courant–Friedrichs–Lewy condition */
	virtual real calculateTimeStep() const override;
	
	/// for tests
	const Mesh* getMesh() const { return mesh; }
	
protected:
	real CourantNumber = 0; ///< number from Courant–Friedrichs–Lewy condition
	GcmMethod gridCharacteristicMethod;
	Corrector* corrector = nullptr;
	InternalOde* internalOde = nullptr;
	SpecialBorder* specialBorder = nullptr;

	Mesh* mesh = nullptr;

	void internalOdeNextStep(const real timeStep);
	void applyCorrectors();
	void moveMesh(const real timeStep);

	USE_AND_INIT_LOGGER("gcm.DefaultSolver")
};



template<class TMesh>
DefaultSolver<TMesh>::
DefaultSolver(const Task& task, AbstractGlobalScene* abstractGlobalScene,
		const GridId gridId_) : Solver(task) {
	GlobalScene* globalScene = dynamic_cast<GlobalScene*>(abstractGlobalScene);
	assert_true(globalScene);
	mesh = new Mesh(task, globalScene, gridId_);
	specialBorder = new SpecialBorder();
}

template<class TMesh>
DefaultSolver<TMesh>::~DefaultSolver() {
	delete specialBorder;
	delete mesh;
}

template<class TMesh>
void DefaultSolver<TMesh>::
beforeStatement(const Statement& statement) {
	CourantNumber = statement.globalSettings.CourantNumber;

	corrector = new Corrector(statement);
	internalOde = new InternalOde(statement);

	specialBorder->beforeStatement(statement);
	mesh->beforeStatement(statement);
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
beforeStage() {
	gridCharacteristicMethod.beforeStage(*mesh);
}

template<class TMesh>
void DefaultSolver<TMesh>::
contactStage(const int s, const real timeStep) {
	gridCharacteristicMethod.contactStage(s, timeStep, *mesh);
}

template<class TMesh>
void DefaultSolver<TMesh>::
privateStage(const int s, const real timeStep) {
	LOG_DEBUG("Start stage " << s << " ... ");
	DATA_BUS::exchangeNodesWithNeighbors(mesh);
	
	specialBorder->applyBorderBeforeStage(mesh, timeStep, s);
	gridCharacteristicMethod.stage(s, timeStep, *mesh);
			//< now actual PDE values is in pdeVariablesNew
	specialBorder->applyBorderAfterStage(mesh, timeStep, s);
	
	// return actual PDE values back to pdeVariables
	std::swap(mesh->pdeVariables, mesh->pdeVariablesNew);
}

template<class TMesh>
void DefaultSolver<TMesh>::
afterStages(const real timeStep) {
	internalOdeNextStep(timeStep);
	moveMesh(timeStep);
	applyCorrectors();
}


template<class TMesh>
void DefaultSolver<TMesh>::
internalOdeNextStep(const real timeStep) {
	if (InternalOde::NonTrivial) {
		assert_eq(mesh->pdeVariables.size(), mesh->odeVariables.size());
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
