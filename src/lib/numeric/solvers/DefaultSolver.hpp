#ifndef LIBGCM_DEFAULTSOLVER_HPP
#define LIBGCM_DEFAULTSOLVER_HPP

#include <lib/numeric/solvers/Solver.hpp>
#include <lib/numeric/gcmethod/GridCharacteristicMethod.hpp>
#include <lib/numeric/border_conditions/BorderConditions.hpp>
#include <lib/util/Logging.hpp>
#include <lib/util/task/Task.hpp>
#include <lib/mesh/AbstractGrid.hpp>
#include <lib/mesh/DataBus.hpp>
#include <lib/rheology/correctors/correctors.hpp>

namespace gcm {
	/**
	 * Class for handling complete time step
	 */
	template<class TMesh>
	class DefaultSolver : public Solver {
	public:
		typedef typename TMesh::Model    Model;
		typedef typename TMesh::Grid     Grid;		
		
		virtual void initializeImpl(const Task& task) override;
		virtual void nextTimeStepImpl() override;
		~DefaultSolver();

		virtual real calculateTau() const override;

		virtual AbstractGrid* getGrid() const {
			return mesh;
		};

	protected:
		real CourantNumber = 0.0; // number from Courant–Friedrichs–Lewy condition
		bool splittingSecondOrder = false; // use or not time second order approach in splitting method

		GridCharacteristicMethod<TMesh>* method = nullptr;
		typename Model::Corrector* corrector = nullptr;
		typename Model::InternalOde* internalOde = nullptr;
		BorderConditions<Model, Grid> borderConditions;

		TMesh* mesh = nullptr;

		USE_AND_INIT_LOGGER("gcm.DefaultSolver");

		void stage(const int s, const real timeStep);
		void internalOdeNextStep(const real timeStep);
		void applyCorrectors();

	};
}

#endif // LIBGCM_DEFAULTSOLVER_HPP
