#ifndef LIBGCM_TASK_HPP
#define LIBGCM_TASK_HPP

#include <map>
#include <memory>
#include <vector>

#include <lib/util/Types.hpp>
#include <lib/util/Concepts.hpp>
#include <lib/util/areas/SphereArea.hpp>
#include <lib/rheology/materials/IsotropicMaterial.hpp>
#include <lib/rheology/materials/OrthotropicMaterial.hpp>

namespace gcm {
	/**
	 * Properties, conditions, tasks in format of the program.
	 * The aim is to separate parsing xml or whatever from program initialization.
	 * Any parser just creates an object of this class and return it.
	 */
	class Task {

	public:
		real xLength = 0.0;
		real yLength = 0.0;
		real zLength = 0.0;

		int accuracyOrder = 0; // order of accuracy of spatial interpolation

		int X = 1; // number of nodes along x direction
		int Y = 1; // number of nodes along y direction
		int Z = 1; // number of nodes along z direction

		real startX = 0.0; // global x-coordinate of the first real node of the grid
		real startY = 0.0; // global y-coordinate of the first real node of the grid
		real startZ = 0.0; // global z-coordinate of the first real node of the grid

		IsotropicMaterial material;

		real CourantNumber = 0.0; // number from Courant–Friedrichs–Lewy condition

		int numberOfSnaps = 0; // how many snaps to calculate
		real T = 0.0; // optional, required time if (numberOfSnaps == 0)

		bool splittingSecondOrder = false;
		bool enableSnapshotting = false;
		bool forceSequence = false; // if true make the grid thinking that the number of
		// processes is one, even if it isn't so actually (for testing purposes)

		/**
		 * All listed here initial conditions will be applied in sequence, not rewriting but adding to each other
		 */
		struct InitialCondition {

			// initial conditions in terms of vector of node values
			struct Vector {
				std::shared_ptr<Area> area;
				std::initializer_list<real> list;
			};
			std::vector<Vector> vectors = {};

			// initial conditions in terms of waves
			struct Wave {
				std::shared_ptr<Area> area;
				Waves::T waveType;
				int direction;
				PhysicalQuantities::T quantity; // quantity to calibrate wave amplitude
				real quantityValue; // value of calibration quantity
			};
			std::vector<Wave> waves = {};

			// initial conditions in terms of physical quantities
			struct Quantity {
				std::shared_ptr<Area> area;
				PhysicalQuantities::T physicalQuantity;
				real value;
			};
			std::vector<Quantity> quantities = {};

		} initialCondition;

		// border conditions for cubic body
		std::map<CUBIC_BORDERS, BorderCondition::T> borderConditions = {
				{CUBIC_BORDERS::X_LEFT,  BorderCondition::T::NON_REFLECTION},
				{CUBIC_BORDERS::X_RIGHT, BorderCondition::T::NON_REFLECTION},
				{CUBIC_BORDERS::Y_LEFT,  BorderCondition::T::NON_REFLECTION},
				{CUBIC_BORDERS::Y_RIGHT, BorderCondition::T::NON_REFLECTION},
				{CUBIC_BORDERS::Z_LEFT,  BorderCondition::T::NON_REFLECTION},
				{CUBIC_BORDERS::Z_RIGHT, BorderCondition::T::NON_REFLECTION}
		};

		Task();
	};
}

#endif // LIBGCM_TASK_HPP
