#ifndef LIBGCM_DEFAULTSOLVER_HPP
#define LIBGCM_DEFAULTSOLVER_HPP

#include <lib/numeric/solvers/Solver.hpp>
#include <lib/numeric/gcmethod/GridCharacteristicMethod.hpp>
#include <lib/util/Logging.hpp>
#include <lib/util/task/Task.hpp>
#include <lib/grid/Grid.hpp>
#include <lib/rheology/correctors/IdealPlasticFlowCorrector.hpp>

namespace gcm {
	/**
	 * Class for handling complete time step
	 */
	template<class TGrid>
	class DefaultSolver : public Solver {
	public:
		virtual void initializeImpl(const Task& task) override;
		virtual void nextTimeStepImpl() override;
		~DefaultSolver();

		virtual real calculateTau() const override;

		virtual Grid* getGrid() const {
			return mesh;
		};

	protected:
		real CourantNumber = 0.0; // number from Courant–Friedrichs–Lewy condition
		bool splittingSecondOrder = false; // use or not time second order approach in splitting method

		GridCharacteristicMethod<TGrid>* method = nullptr;
		IdealPlasticFlowCorrector<TGrid>* plasticFlowCorrector = nullptr;

		/** After calculation of step, actual values are always here */
		TGrid* mesh = nullptr;
		/** This is auxiliary mesh */
		TGrid* newMesh = nullptr;

		USE_AND_INIT_LOGGER("gcm.DefaultSolver");

		void stage(const int s, const real timeStep);

	};
}

#endif // LIBGCM_DEFAULTSOLVER_HPP
