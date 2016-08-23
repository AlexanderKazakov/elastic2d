#ifndef LIBGCM_TASK_HPP
#define LIBGCM_TASK_HPP


#include <map>
#include <memory>
#include <vector>

#include <lib/util/Types.hpp>
#include <lib/util/Enum.hpp>
#include <lib/util/Area.hpp>
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
	
	std::string id; ///< name of statement
	
	
	struct GlobalSettings {
		real CourantNumber = 0; ///< number in Courant–Friedrichs–Lewy condition
		int numberOfSnaps = 0;
		int stepsPerSnap = 1;
		real requiredTime = 0;  ///< optional, required time if (numberOfSnaps <= 0)
	} globalSettings;
	
	
	struct MaterialCondition {
		typedef std::shared_ptr<AbstractMaterial> Material;
		
		enum class Type {
			BY_AREAS,
			BY_BODIES,
		} type = Type::BY_AREAS;
		
		
		struct ByAreas {
			Material defaultMaterial; ///< everywhere but inhomogenities
			struct Inhomogenity {
				std::shared_ptr<Area> area;
				Material material;
			};
			std::vector<Inhomogenity> materials; ///< inhomogenities
		} byAreas;
		
		
		struct ByBodies {
			typedef size_t BodyId;
			std::map<BodyId, Material> bodyMaterialMap;
		} byBodies;
		
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
	
	
	struct CubicGridBorderCondition {
		std::shared_ptr<Area> area;
		std::map<PhysicalQuantities::T, TimeDependency> values;
	};
	std::vector<CubicGridBorderCondition> cubicGridBorderConditions;
	
	
	struct BorderCondition {
		std::shared_ptr<Area> area;
		BorderConditions::T type;
		std::vector<TimeDependency> values; ///< size must be equal to 
		///< number of outer characteristics
	};
	std::vector<BorderCondition> borderConditions;
	
	
	struct Fracture {
		int direction; ///< crossing axis
		int index; ///< index at crossing axis
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
};


/**
 * Struct to initialize the program.
 * 
 * Note: However almost classes initialized by sending to them (const Task& task),
 * they must use for initialization only its own components of Task,
 * because other components may not be initialized.
 * For example, some special snapshotter for CubicGrid must not use task.cubicGrid,
 * because from cubicGrid.lenghts, cubicGrid.h, cubicGrid.sizes only two 
 * can be initialized at time.
 */
struct Task {
	
	struct Body {
		Materials::T materialId; ///< body material
		Models::T modelId;       ///< body rheology
	};
	/// Bodies sorted by unique id
	std::map<size_t, Body> bodies;
	
	
	struct GlobalSettings {
		/// Space dimensionality
		int dimensionality;
		/// Type of grid used in calculations
		Grids::T gridId;
		/// If true, independently from number of working processes
		/// calculation of a concrete statement is performed in sequence.
		/// Useful for inverse problem calculation and MPI testing.
		bool forceSequence = false;
		/// On/off deformations and bodies motion
		bool movable = false;
		/// List of snapshotters to use
		std::vector<Snapshotters::T> snapshottersId;		
	} globalSettings;
	
	
	struct CubicGrid {
		std::vector<real> lengths;  ///< lengthes of cube in each direction
		std::vector<real> h;        ///< spatial steps in each direction
		std::vector<int>  sizes;    ///< number of nodes along each direction
		std::vector<real> startR;   ///< global coordinates of the first real node
		int borderSize = 0;         ///< number of virtual border nodes for one border point
		bool forceSequence = false; ///< behave such as Mpi::Size() == 1
	} cubicGrid;
	
	
	struct SimplexGrid {
		
		enum class Mesher {
			CGAL_MESHER,
			INM_MESHER,
		} mesher = Mesher::CGAL_MESHER;
		
		/// effective spatial step for mesher
		real spatialStep = 0;
		
		/// option for Cgal3DMesher only - use true for figures with sharp edges
		bool detectSharpEdges = false; 
		
		/// file with some initial data for mesher
		std::string fileName;
		
		/// denominator to scale the points after meshing
		real scale = 1;
		
		/// for Cgal2DMesher only
		struct Body {
			typedef std::array<real, 2> Point;
			typedef std::vector<Point> Border;
			
			size_t id;                  ///< body indicator > 0 @see Task::Body
			Border outer;               ///< outer border of the body
			std::vector<Border> inner;  ///< borders of the inner cavities
		};
		std::vector<Body> bodies; ///< list of bodies contained in the grid
	
	} simplexGrid;
	
	
	struct ContactCondition {
		typedef std::pair<size_t, size_t> GridsContact;
		ContactConditions::T defaultCondition;
		std::map<GridsContact, ContactConditions::T> gridToGridConditions;
	} contactCondition;
	
	
	/// list of statements to calculate on the same geometry
	std::vector<Statement> statements;
	
};

/** @} */
}

#endif // LIBGCM_TASK_HPP
