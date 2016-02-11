#ifndef LIBGCM_TEST_WRAPPERS_HPP
#define LIBGCM_TEST_WRAPPERS_HPP

#include <lib/Engine.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/grid/StructuredGrid.hpp>

namespace gcm {

	/**
	 * Wrappers for access to protected members for testing purposes
	 */

	template<class TModel>
	class StructuredGridWrapper : public StructuredGrid<TModel> {
	public:
		linal::Vector<3> getHForTest() const { return this->h; };
		linal::VectorInt<3> getSizesForTest() const { return this->sizes; };
		int getStartXForTest() const { return this->globalStartXindex; };

		typename StructuredGrid<TModel>::Matrix interpolateValuesAroundForTest
				(const int stage, const int x, const int y, const int z,
				 const typename StructuredGrid<TModel>::PdeVector& dx) const {
			return this->interpolateValuesAround(stage, x, y, z, dx);
		};

		void findSourcesForInterpolationForTest
				(const int stage, const int x, const int y, const int z, const real &dx,
				 std::vector<typename StructuredGrid<TModel>::PdeVector>& src) const {
			return this->findSourcesForInterpolation(stage, x, y, z, dx, src);
		};

		void changeRheology(const real &rho2rho0, const real &lambda2lambda0, const real &mu2mu0) {
			IsotropicMaterial oldMaterial = this->getMatrix(0, 0, 0)->getMaterial();
			IsotropicMaterial newMaterial
					(rho2rho0 * oldMaterial.rho, lambda2lambda0 * oldMaterial.lambda, mu2mu0 * oldMaterial.mu);
			auto newRheologyMatrix = new typename TModel::GCM_MATRICES(newMaterial);

			for (int x = 0; x < this->sizes(0); x++) {
				for (int y = 0; y < this->sizes(1); y++) {
					for (int z = 0; z < this->sizes(2); z++) {
						if (y * this->h(1) >= 0.5) {
							this->gcmMatrices[this->getIndex(x, y, z)] = newRheologyMatrix;
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
			return static_cast<StructuredGridWrapper<typename TGrid::Model>*>(this->grid);
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
