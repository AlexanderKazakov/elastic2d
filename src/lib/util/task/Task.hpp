#ifndef LIBGCM_TASK_HPP
#define LIBGCM_TASK_HPP

#include <map>
#include <memory>
#include <vector>

#include <lib/linal/special/VectorInt.hpp>
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
		int accuracyOrder = 0; // order of accuracy of spatial interpolation

		linal::Vector<3> lengthes = {0, 0, 0}; // lengthes of cube in each direction
		linal::VectorInt<3> sizes = {1, 1, 1}; // number of nodes along each direction
		linal::Vector<3> startR = {0, 0, 0}; // global coordinates of the first real node

		IsotropicMaterial material;

		real CourantNumber = 0.0; // number from Courant–Friedrichs–Lewy condition
		int numberOfSnaps = 0; // how many snaps to calculate
		real T = 0.0; // optional, required time if (numberOfSnaps == 0)

		bool splittingSecondOrder = false;
		bool enableSnapshotting = false;
		bool forceSequence = false; // if true make the grid thinking that the number of
		// processes is one, even if it isn't so actually

		bool plasticityFlowCorrector = false;

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
	};
}

#endif // LIBGCM_TASK_HPP
