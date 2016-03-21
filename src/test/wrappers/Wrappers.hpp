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
		typedef typename TMesh::Model Model;
		typedef typename TMesh::GcmMatricesPtr GcmMatricesPtr;
		typedef typename TMesh::GCM_MATRICES GCM_MATRICES;
		typedef typename TMesh::PdeVector PdeVector;
		typedef typename TMesh::Iterator Iterator;

		MeshWrapper(const Task& task) : TMesh(task) { }
	};

	template<class TGrid>
	class DefaultSolverWrapper : public DefaultSolver<TGrid> {
	public:
		void stageForTest(const int s, const real timeStep) { return this->stage(s, timeStep); }
		real getTauForTest() const { return this->calculateTimeStep(); }

		MeshWrapper<TGrid>* getMesh() const {
			return static_cast<MeshWrapper<TGrid>*>(this->mesh);
		}
	};

	template<class TGrid>
	class EngineWrapper : public Engine {
	public:
		DefaultSolverWrapper<TGrid>* getSolverForTest() const {
			return static_cast<DefaultSolverWrapper<TGrid>*>(this->solver);
		}
		void runStatementForTest() { this->runStatement(); }
		
		void beforeStatementForTest(const Statement& statement) {
			this->beforeStatement(statement);
		}
	};
}

#endif // LIBGCM_TEST_WRAPPERS_HPP
