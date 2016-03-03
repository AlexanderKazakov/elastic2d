#ifndef LIBGCM_EXTERNALTASK_HPP
#define LIBGCM_EXTERNALTASK_HPP

#include <map>
#include <memory>
#include <vector>

#include <lib/linal/linal.hpp>
#include <lib/util/Types.hpp>
#include <lib/util/Enum.hpp>
#include <lib/util/areas/areas.hpp>
#include <lib/rheology/materials/materials.hpp>

namespace gcm {
	/**
	 * Properties, conditions, tasks in external format - strings and numbers.
	 * The aim is to separate parsing xml or whatever from program initialization.
	 * Any parser just creates an object of this class and returns it.
	 * Then it's translated to Task - program's format of conditions.
	 */
	
	/** 
	 * Set of all conditions except mesh geometry.
	 * This is part of task for one statement.
	 */
	struct Statement {
		typedef std::function<real(real)> TimeDependency;

		std::string id; // name of statement

		real CourantNumber = 0.0; // number from Courant–Friedrichs–Lewy condition
		int numberOfSnaps = 0;
		int stepsPerSnap = 1;
		real T = 0.0; // optional, required time if (numberOfSnaps == 0)

		IsotropicMaterial isotropicMaterial;
		OrthotropicMaterial orthotropicMaterial;

		/**
		 * All listed here initial conditions will be applied in sequence, 
		 * not rewriting but adding to each other
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

		struct BorderCondition {
			std::shared_ptr<Area> area;
			std::map<PhysicalQuantities::T, TimeDependency> values;
		};
		std::vector<BorderCondition> borderConditions = {};

		struct Fracture {
			int direction;
			real coordinate;
			std::shared_ptr<Area> area;
			std::map<PhysicalQuantities::T, TimeDependency> values;
		};
		std::vector<Fracture> fractures = {};

		std::vector<PhysicalQuantities::T> quantitiesToVtk = {};

		struct Detector {
			std::vector<PhysicalQuantities::T> quantities = {};
			std::shared_ptr<Area> area;
		} detector;
	};
	
	struct Task {
		linal::Vector<3> lengthes = {0, 0, 0}; // lengthes of cube in each direction
		linal::VectorInt<3> sizes = {1, 1, 1}; // number of nodes along each direction
		linal::Vector<3> startR = {0, 0, 0}; // global coordinates of the first real node
		real spatialStep = 0; // effective spatial step for unstructured grids

		int accuracyOrder = 0; // order of accuracy of spatial interpolation
		bool splittingSecondOrder = false;
		bool enableSnapshotting = false;
		bool forceSequence = false; // if true make meshes thinking that the number of
		// processes is one, even if it isn't so actually

		// list of statements to calculate on the same geometry
		std::vector<Statement> statements;
	};
}

#endif // LIBGCM_EXTERNALTASK_HPP
