#ifndef LIBGCM_TASK_HPP
#define LIBGCM_TASK_HPP

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
 * @defgroup task Task and statements
 * Properties, conditions, tasks in native format of the program.
 * Used in initialization of the program.
 * @{
 */

/**
 * Set of all conditions for one statement.
 * Task can contain several statements.
 */
struct Statement {
	typedef std::function<real(real)> TimeDependency;

	std::string id;         ///< name of statement

	struct GlobalSettings {
		real CourantNumber = 0.0;        ///< number from Courant–Friedrichs–Lewy condition
		int numberOfSnaps = 0;
		int stepsPerSnap = 1;
		real requiredTime = 0.0;         ///< optional, required time if (numberOfSnaps <= 0)
	} globalSettings;

	struct MaterialCondition {
		std::shared_ptr<AbstractMaterial> defaultMaterial;
		struct Inhomogenity {
			std::shared_ptr<Area> area;
			std::shared_ptr<AbstractMaterial> material;
		};
		std::vector<Inhomogenity> materials;
	} materialConditions;

	/**
	 * All listed here initial conditions will be applied in sequence,
	 * not rewriting but adding to each other
	 */
	struct InitialCondition {
		/// initial conditions in terms of vector of node values
		struct Vector {
			std::shared_ptr<Area> area;
			std::initializer_list<real> list;
		};
		std::vector<Vector> vectors;

		/// initial conditions in terms of waves
		struct Wave {
			std::shared_ptr<Area> area;
			Waves::T waveType;
			int direction;
			PhysicalQuantities::T quantity; ///< quantity to calibrate wave amplitude
			real quantityValue;             ///< value of calibration quantity
		};
		std::vector<Wave> waves;

		/// initial conditions in terms of physical quantities
		struct Quantity {
			std::shared_ptr<Area> area;
			PhysicalQuantities::T physicalQuantity;
			real value;
		};
		std::vector<Quantity> quantities;

	} initialCondition;

	struct BorderCondition {
		std::shared_ptr<Area> area;
		std::map<PhysicalQuantities::T, TimeDependency> values;
	};
	std::vector<BorderCondition> borderConditions;

	struct Fracture {
		int direction;
		real coordinate;
		std::shared_ptr<Area> area;
		std::map<PhysicalQuantities::T, TimeDependency> values;
	};
	std::vector<Fracture> fractures;

	struct VtkSnapshotter {
		/// in order to not dump big vtk snaps every statement
		/// but just sometimes for eye-checking
		bool enableSnapshotting = false;
		/// list of physical quantities to write to vtk
		std::vector<PhysicalQuantities::T> quantitiesToSnap;
	} vtkSnapshotter;

	struct Detector {
		std::vector<PhysicalQuantities::T> quantities;
		std::shared_ptr<Area> area;
	} detector;

	struct Binary2DSeismograph {
		PhysicalQuantities::T quantityToWrite;
	} binary2DSeismograph;
};

/**
 * Struct to initialize the program.
 */
struct Task {

	Materials::T materialId;
	Grids::T gridId;
	Models::T modelId;
	std::vector<Snapshotters::T> snapshottersId;

	struct GlobalSettings {
		/// If true, independently from number of working processes
		/// calculation of a concrete statement is performed in sequence.
		/// Useful for inverse problem calculation and MPI testing.
		bool forceSequence = false;
	} globalSettings;

	struct CubicGrid {
		int dimensionality = 0;     ///< spatial dimensionality of the grid 
		///< (model can have different)
		Real3 lengths = {0, 0, 0};  ///< lengthes of cube in each direction
		Real3 h = {0, 0, 0};        ///< spatial steps in each direction
		Int3 sizes = {0, 0, 0};     ///< number of nodes along each direction
		Real3 startR = {0, 0, 0};   ///< global coordinates of the first real node
		int borderSize = 0;         ///< number of virtual border nodes for one border point
	} cubicGrid;

	struct Cgal2DGrid {
		bool movable = false;       ///< deformable(true) or immutable(false) grid
		real spatialStep = 0;       ///< effective spatial step
		
		struct Body {
			typedef std::vector<Real2> Border;
			Body(const Border& outer_, const std::vector<Border>& inner_) :
				outer(outer_), inner(inner_) { }
			
			Border outer;               ///< outer border of the body
			std::vector<Border> inner;  ///< borders of the inner cavities
		};
		std::vector<Body> bodies; ///< list of bodies contained in the grid
		
	} cgal2DGrid;

	/// list of statements to calculate on the same geometry
	std::vector<Statement> statements;
};

/** @} */
}

#endif // LIBGCM_TASK_HPP
