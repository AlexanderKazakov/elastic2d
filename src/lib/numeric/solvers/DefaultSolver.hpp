#ifndef LIBGCM_DEFAULTSOLVER_HPP
#define LIBGCM_DEFAULTSOLVER_HPP

#include <lib/numeric/solvers/Solver.hpp>
#include <lib/numeric/gcm/GridCharacteristicMethod.hpp>
#include <lib/numeric/border_conditions/BorderConditions.hpp>
#include <lib/util/Logging.hpp>
#include <lib/util/task/Task.hpp>
#include <lib/mesh/grid/AbstractGrid.hpp>
#include <lib/mesh/DataBus.hpp>
#include <lib/mesh/MeshMover.hpp>
#include <lib/rheology/correctors/correctors.hpp>

namespace gcm {
	/**
	 * Class for handling complete time step
	 */
	template<class TMesh>
	class DefaultSolver : public Solver {
	public:
		typedef typename TMesh::Model                          Model;
		typedef typename TMesh::Grid                           Grid;
		typedef typename TMesh::Material                       Material;
		typedef typename Model::Corrector                      Corrector;
		typedef typename Model::InternalOde                    InternalOde;
		typedef BorderConditions<Model, Grid, Material>        Border;
		typedef DataBus<Model, Grid, Material>                 DATA_BUS;
		typedef MeshMover<Model, Grid, Material>               MESH_MOVER;
		typedef GridCharacteristicMethod<TMesh>                GCM;

		virtual ~DefaultSolver();
		virtual real calculateTau() const override;
		virtual AbstractGrid* getGrid() const { return mesh; }

	protected:
		real CourantNumber = 0.0; // number from Courant–Friedrichs–Lewy condition
		bool splittingSecondOrder = false; // time second order approach in splitting method

		Corrector* corrector = nullptr;
		InternalOde* internalOde = nullptr;
		Border* borderConditions = nullptr;

		TMesh* mesh = nullptr;
		
		virtual void initializeImpl(const Task& task) override;
		virtual void beforeStatementImpl(const Statement& statement) override;
		virtual void nextTimeStepImpl() override;
		virtual void afterStatementImpl() override;

		void stage(const int s, const real timeStep);
	private:
		void internalOdeNextStep(const real timeStep);
		void applyCorrectors();
		void moveMesh(const real timeStep);

		USE_AND_INIT_LOGGER("gcm.DefaultSolver")
	};
}

#endif // LIBGCM_DEFAULTSOLVER_HPP
