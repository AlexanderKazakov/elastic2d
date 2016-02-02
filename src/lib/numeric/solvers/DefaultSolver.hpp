#ifndef LIBGCM_DEFAULTSOLVER_HPP
#define LIBGCM_DEFAULTSOLVER_HPP

#include <lib/numeric/solvers/Solver.hpp>
#include <lib/numeric/gcmethod/GridCharacteristicMethod.hpp>
#include <lib/util/Logging.hpp>
#include <lib/util/task/Task.hpp>
#include <lib/grid/Grid.hpp>
#include <lib/rheology/correctors/correctors.hpp>

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
		/*IdealPlasticFlowCorrector<TGrid>*/
		typename TGrid::Model::Corrector* corrector = nullptr;
		typename TGrid::Model::InternalOde* internalOde = nullptr;

		/** Actual values are always here */
		TGrid* mesh = nullptr;
		/** This is auxiliary mesh */
		TGrid* newMesh = nullptr;
		bool odeShiftedFromPde = false;

		USE_AND_INIT_LOGGER("gcm.DefaultSolver");

		void stage(const int s, const real timeStep);
		/**
		 * Make the most actual at the moment PDE values
		 * to be at the same mesh with the most actual at the moment ODE values.
		 * TODO - remove this brainfuck
		 * @warning this function fixes situation just for now,
		 * when some more complicated splitting or rheology or coordinates appears it will be broken
		 */
		void fixVariablesOrder();
		void internalOdeNextStep(const real timeStep);
		void applyCorrectors();

	};
}

#endif // LIBGCM_DEFAULTSOLVER_HPP
