#ifndef LIBGCM_TASK_HPP
#define LIBGCM_TASK_HPP


#include <map>
#include <memory>
#include <vector>

#include <libgcm/util/infrastructure/infrastructure.hpp>
#include <libgcm/util/Enum.hpp>
#include <libgcm/util/math/Area.hpp>
#include <libgcm/rheology/materials/materials.hpp>


namespace gcm {

/**
 * Struct to initialize the program.
 * Properties, conditions, tasks in native format.
 * @note However almost classes initialized by sending to them (const Task& task),
 * they must use for initialization only its own components of Task,
 * because other components may not be initialized or valid.
 */
struct Task {
	typedef std::function<real(real)> TimeDependency;
	
	struct Body {
		/// Body material
		Materials::T materialId;
		/// Body rheology
		Models::T modelId;
		/// List of ODE, correctors and other non-wave features to apply
		std::vector<Odes::T> odes;
	};
	/// Bodies sorted by unique id
	std::map<size_t, Body> bodies;
	
	
	struct GlobalSettings {
		/// Space dimensionality
		int dimensionality;
		/// Type of grid used in calculations
		Grids::T gridId;
		/// If true, independently from number of working processes
		/// calculation is performed in sequence. Useful for MPI testing.
		bool forceSequence = false;
		/// List of snapshotters to use
		std::vector<Snapshotters::T> snapshottersId;
		/// Name of the folder to write snapshots into relative to ./snapshots/
		std::string outputDirectory = "";
		/// Number in Courant–Friedrichs–Lewy condition
		real CourantNumber = 0;
		/// Number of snapshots to make in calculations
		int numberOfSnaps = 0;
		/// Total number of time steps == numberOfSnaps * stepsPerSnap
		int stepsPerSnap = 1;
		/// Optional, required time if (numberOfSnaps <= 0)
		real requiredTime = 0;
		/// On/off logging from engine at each time step
		bool verboseTimeSteps = true;
		/// Type of gcm-method to use (only for simplex engine now)
		GcmType gcmType = GcmType::ADVECT_RIEMANN_INVARIANTS;
		/// Type of splitting by directions (only for simplex engine now)
		SplittingType splittingType = SplittingType::PRODUCT;
	} globalSettings;
	
	
	struct CubicGrid {
		struct Cube {
			/// number of nodes along each direction
			std::vector<int> sizes;
			/// global index of the most left real node
			std::vector<int> start;
		};
		/// spatial steps in each coordinate direction
		std::vector<real> h;
		/// number of ghost border nodes used for border and contact calculation
		int borderSize;
		/// list of cubic bodies sorted by unique id @see Task::Body
		std::map<size_t, Cube> cubics;
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
		
		/// On/off deformations and bodies motion
		bool movable = false;
		
		/// Method of border and contact nodes calculation
		BorderCalcMode borderCalcMode = BorderCalcMode::GLOBAL_BASIS;
		
		
		/// for Cgal2DMesher only @{
		struct Body {
			typedef std::array<real, 2> Point;
			typedef std::vector<Point> Border;
			
			size_t id;                  ///< body indicator > 0 @see Task::Body
			Border outer;               ///< outer border of the body
			std::vector<Border> inner;  ///< borders of the inner cavities
		};
		std::vector<Body> bodies; ///< list of bodies contained in the grid
		/// @}
		
	} simplexGrid;
	
	
	/// Calculation basis to use (for simplex engine only).
	/// If not specified, new random basis is used at every time step.
	/// The components of matrix are stored in C-style (string-by-string).
	/// The direction of calculation on i'th stage is i'th column of the matrix.
	std::vector<real> calculationBasis;
	
	
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
	
	
	struct CubicBorderCondition {
		int direction;
		std::shared_ptr<Area> area;
		typedef std::map<PhysicalQuantities::T, TimeDependency> Values;
		Values values;
	};
	/// lists of conditions for all bodies;
	/// conditions in lists are applied in sequence overwriting each other,
	/// i.e if two conditions share some node, the last one will be applied
	std::map<size_t, std::vector<CubicBorderCondition>> cubicBorderConditions;
	
	
	struct BorderCondition {
		std::shared_ptr<Area> area;
		/// we calculate multicontact nodes as borders, but sometimes
		/// it can be not appropriate to apply some border conditions to them
		bool useForMulticontactNodes = true;
		BorderConditions::T type;
		///< size must be equal to number of outer characteristics
		std::vector<TimeDependency> values;
	};
	std::vector<BorderCondition> borderConditions;
	
	
	struct ContactCondition {
		typedef std::pair<size_t, size_t> GridsContact;
		ContactConditions::T defaultCondition;
		std::map<GridsContact, ContactConditions::T> gridToGridConditions;
	} contactCondition;
	
	
	struct VtkSnapshotter {
		/// list of physical quantities to write to vtk
		std::vector<PhysicalQuantities::T> quantitiesToSnap;
	} vtkSnapshotter;
	
	struct Detector {
		std::vector<PhysicalQuantities::T> quantities;
		std::shared_ptr<Area> area;
		size_t gridId;
	} detector;
	
};


}

#endif // LIBGCM_TASK_HPP
