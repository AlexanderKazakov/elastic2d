#ifndef LIBGCM_TEST_WRAPPERS_HPP
#define LIBGCM_TEST_WRAPPERS_HPP

#include <lib/Engine.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>

namespace gcm {
/**
 * Wrappers for access to protected members for testing purposes
 */

template<class TMesh>
class MeshWrapper : public TMesh {
public:
	typedef typename TMesh::Model          Model;
	typedef typename TMesh::GcmMatricesPtr GcmMatricesPtr;
	typedef typename TMesh::GCM_MATRICES   GCM_MATRICES;
	typedef typename TMesh::PdeVector      PdeVector;
	typedef typename TMesh::Iterator       Iterator;

	MeshWrapper(const Task& task) : TMesh(task) { }

	void beforeStatementForTest(const Statement& statement) {
		this->beforeStatement(statement);
	}

};


template<class TMesh>
class DefaultSolverWrapper : public DefaultSolver<TMesh> {
public:
	DefaultSolverWrapper(const Task& task) : DefaultSolver<TMesh>(task) { }

	// TODO - move to GcmMethod test
	void stageForTest(const int s, const real timeStep) {
		this->stage(s, timeStep);
	}

	MeshWrapper<TMesh>* getMesh() {
		return static_cast<MeshWrapper<TMesh>*>(this->getActualMesh());
	}

};


template<class TMesh>
class EngineWrapper : public Engine {
public:
	EngineWrapper(const Task& task_) : Engine(task_) { }

	DefaultSolverWrapper<TMesh>* getSolverForTest() const {
		return static_cast<DefaultSolverWrapper<TMesh>*>(this->solver);
	}

	void runStatementForTest() {
		this->runStatement();
	}

	void beforeStatementForTest(const Statement& statement) {
		this->beforeStatement(statement);
	}

};


}

#endif // LIBGCM_TEST_WRAPPERS_HPP
