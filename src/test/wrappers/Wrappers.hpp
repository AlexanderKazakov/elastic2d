#ifndef LIBGCM_TEST_WRAPPERS_HPP
#define LIBGCM_TEST_WRAPPERS_HPP

#include <lib/Engine.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/grid/StructuredGrid.hpp>

namespace gcm {
	/**
	 * Wrappers for access to protected members for testing purposes
	 */

	template<class TMesh>
	class MeshWrapper : public TMesh {
	public:
		typedef typename TMesh::GCM_MATRICES GCM_MATRICES;
		typedef typename TMesh::PdeVector PdeVector;
		typedef typename TMesh::Iterator Iterator;

		PdeVector getPde(const int x, const int y, const int z) const {
			auto it = Iterator({x, y, z}, this->sizes);
			return this->pde(it);
		};
		GCM_MATRICES* getMatrix(const int x, const int y, const int z) const {
			return this->matrix(Iterator({x, y, z}));
		};
		GCM_MATRICES*& getMatrix(const int x, const int y, const int z) {
			return this->_matrix(Iterator({x, y, z}));
		};

		void changeRheology(const real rho2rho0, const real lambda2lambda0, const real mu2mu0) {
			IsotropicMaterial oldMaterial = this->getMatrix(0, 0, 0)->getMaterial();
			IsotropicMaterial newMaterial
					(rho2rho0 * oldMaterial.rho, lambda2lambda0 * oldMaterial.lambda, mu2mu0 * oldMaterial.mu);
			auto newRheologyMatrix = new GCM_MATRICES(newMaterial);

			for (int x = 0; x < this->sizes(0); x++) {
				for (int y = 0; y < this->sizes(1); y++) {
					for (int z = 0; z < this->sizes(2); z++) {
						if (y * this->h(1) >= 0.5) {
							this->getMatrix(x, y, z) = newRheologyMatrix;
						}
					}
				}
			}

			this->maximalLambda = fmax(this->maximalLambda, newRheologyMatrix->getMaximalEigenvalue());
		}
	};

	template<class TGrid>
	class DefaultSolverWrapper : public DefaultSolver<TGrid> {
	public:
		void stageForTest(const int s, const real &timeStep) { return this->stage(s, timeStep); }
		real getTauForTest() const { return this->calculateTau(); };

		MeshWrapper<TGrid>* getMesh() const {
			return static_cast<MeshWrapper<TGrid>*>(this->grid);
		};
	};

	template<class TGrid>
	class EngineWrapper : public Engine {
	public:
		DefaultSolverWrapper<TGrid>* getSolverForTest() const {
			return static_cast<DefaultSolverWrapper<TGrid>*>(this->solver);
		}
	};
}

#endif // LIBGCM_TEST_WRAPPERS_HPP
