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
		typedef TMesh                                         Mesh;
		typedef typename Mesh::Model                          Model;
		typedef typename Mesh::Grid                           Grid;
		typedef typename Mesh::Material                       Material;

		typedef typename Model::Corrector                     Corrector;
		typedef typename Model::InternalOde                   InternalOde;

		typedef BorderConditions<Model, Grid, Material>       Border;
		typedef DataBus<Model, Grid, Material>                DATA_BUS;
		typedef MeshMover<Model, Grid, Material>              MESH_MOVER;
		typedef GridCharacteristicMethod<Mesh>                GCM;

		DefaultSolver(const Task& task);
		virtual ~DefaultSolver();

		virtual void beforeStatement(const Statement& statement) override;
		virtual void afterStatement() override;

		virtual void nextTimeStep(const real timeStep) override;

		/** @return grid with actual values */
		virtual AbstractGrid* getActualGrid() const { return mesh; }

		/** Calculate time step from Courant–Friedrichs–Lewy condition */
		virtual real calculateTimeStep() const override;

	protected:
		real CourantNumber = 0.0; // number from Courant–Friedrichs–Lewy condition
		bool splittingSecondOrder = false; // time second order approach in splitting method

		Corrector* corrector = nullptr;
		InternalOde* internalOde = nullptr;
		Border* borderConditions = nullptr;

		Mesh* mesh = nullptr;

		void stage(const int s, const real timeStep);
		void internalOdeNextStep(const real timeStep);
		void applyCorrectors();
		void moveMesh(const real timeStep);

		USE_AND_INIT_LOGGER("gcm.DefaultSolver")
	};
}

#endif // LIBGCM_DEFAULTSOLVER_HPP
