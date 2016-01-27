#ifndef LIBGCM_TEST_WRAPPERS_HPP
#define LIBGCM_TEST_WRAPPERS_HPP

#include <lib/Engine.hpp>
#include <lib/grid/StructuredGrid.hpp>

namespace gcm {

	/**
	 * Wrappers for access to protected members for testing purposes
	 */

	template<class TModel>
	class StructuredGridWrapper : public StructuredGrid<TModel> {
	public:
		real getH0ForTest() const { return this->h(0); };

		real getH1ForTest() const { return this->h(1); };

		int getYForTest() const { return this->sizes(1); };

		int getXForTest() const { return this->sizes(0); };

		int getStartXForTest() const { return this->globalStartXindex; };

		const typename StructuredGrid<TModel>::Node &getNodeForTest(const int x, const int y, const int z) const {
			return this->get(x, y, z);
		};

		typename StructuredGrid<TModel>::Matrix interpolateValuesAroundForTest
				(const int stage, const int x, const int y, const int z,
				 const typename StructuredGrid<TModel>::Vector& dx) const {
			return this->interpolateValuesAround(stage, x, y, z, dx);
		};

		void findSourcesForInterpolationForTest
				(const int stage, const int x, const int y, const int z, const real &dx,
				 std::vector<typename StructuredGrid<TModel>::Vector>& src) const {
			return this->findSourcesForInterpolation(stage, x, y, z, dx, src);
		};

		void changeRheology(const real &rho2rho0, const real &lambda2lambda0, const real &mu2mu0) {
			IsotropicMaterial oldMaterial = (*this)(0, 0, 0).matrix->getMaterial();
			IsotropicMaterial newMaterial
					(rho2rho0 * oldMaterial.rho, lambda2lambda0 * oldMaterial.lambda, mu2mu0 * oldMaterial.mu);
			auto newRheologyMatrix = std::make_shared<typename TModel::GCM_MATRICES>(newMaterial);

			for (int x = 0; x < this->sizes(0); x++) {
				for (int y = 0; y < this->sizes(1); y++) {
					for (int z = 0; z < this->sizes(2); z++) {
						if (y * this->h(1) >= 0.5) {
							(*this)(x, y, z).matrix = newRheologyMatrix;
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

		StructuredGridWrapper<typename TGrid::Model>* getMesh() const {
			return static_cast<StructuredGridWrapper<typename TGrid::Model>*>(this->mesh);
		};

		StructuredGridWrapper<typename TGrid::Model>* getNewMesh() const {
			return static_cast<StructuredGridWrapper<typename TGrid::Model>*>(this->newMesh);
		};
	};

	template<class TGrid>
	class EngineWrapper : public Engine<TGrid> {
	public:
		DefaultSolverWrapper<TGrid>* getSolver() const {
			return static_cast<DefaultSolverWrapper<TGrid>*>(this->solver);
		}
	};
}

#endif // LIBGCM_TEST_WRAPPERS_HPP
