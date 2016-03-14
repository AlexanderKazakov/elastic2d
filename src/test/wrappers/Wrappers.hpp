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

		void changeRheology(const real rho2rho0, const real lambda2lambda0, const real mu2mu0) {
			IsotropicMaterial oldMaterial = *(this->material({0, 0, 0}));
			IsotropicMaterial newMaterial
					(rho2rho0 * oldMaterial.rho, lambda2lambda0 * oldMaterial.lambda, mu2mu0 * oldMaterial.mu);
			auto newRheologyMatrix = std::make_shared<GCM_MATRICES>();
			Model::constructGcmMatrices(newRheologyMatrix, PdeVector::zeros(), newMaterial);

			for (int x = 0; x < this->sizes(0); x++) {
				for (int y = 0; y < this->sizes(1); y++) {
					if (y * this->h(1) >= 0.5) {
						this->_matrix({x, y, 0}) = newRheologyMatrix;
					}
				}
			}

			this->maximalLambda = fmax(this->maximalLambda, newRheologyMatrix->getMaximalEigenvalue());
		}
	};

	template<class TGrid>
	class DefaultSolverWrapper : public DefaultSolver<TGrid> {
	public:
		void stageForTest(const int s, const real timeStep) { return this->stage(s, timeStep); }
		real getTauForTest() const { return this->calculateTau(); }

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
