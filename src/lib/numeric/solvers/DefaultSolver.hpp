#ifndef LIBGCM_DEFAULTSOLVER_HPP
#define LIBGCM_DEFAULTSOLVER_HPP

#include <lib/numeric/gcmethod/GridCharacteristicMethod.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/util/Logging.hpp>

namespace gcm {
	template <class TNode> class Engine;

	/**
	 * Class for handling complete time step
	 */
	template<class TNode>
	class DefaultSolver {
	public:
		void initialize(const Task& task);
		~DefaultSolver();

		void nextTimeStep();

	protected:
		real currentTime = 0.0;
		real CourantNumber = 0.0; // number from Courant–Friedrichs–Lewy condition
		bool splittingSecondOrder = false; // use or not time second order approach in splitting method

		GridCharacteristicMethod<TNode>* method = nullptr;

		/** After calculation of step, actual values are always here */
		StructuredGrid<TNode>* mesh = nullptr;
		/** This is auxiliary mesh */
		StructuredGrid<TNode>* newMesh = nullptr;

		USE_AND_INIT_LOGGER("gcm.DefaultSolver");

		void stage(const int s, const real &timeStep);

		/** Calculate time step from Courant–Friedrichs–Lewy condition */
		real calculateTau() const;

		friend class Engine<TNode>;
	};
}

#endif // LIBGCM_DEFAULTSOLVER_HPP
