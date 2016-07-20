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
		real CourantNumber = 0; ///< number from Courant–Friedrichs–Lewy condition
		int numberOfSnaps = 0;
		int stepsPerSnap = 1;
		real requiredTime = 0;  ///< optional, required time if (numberOfSnaps <= 0)
	} globalSettings;

	struct MaterialCondition {
		typedef std::shared_ptr<AbstractMaterial> Material;
		
		enum class Type {
			BY_AREAS,
			BY_CELLS,
		} type = Type::BY_AREAS;
		
		/// @name for BY_AREAS @{
		Material defaultMaterial; ///< everywhere but inhomogenities
		struct Inhomogenity {
			std::shared_ptr<Area> area;
			Material material;
		};
		std::vector<Inhomogenity> materials; ///< inhomogenities
		/// @}
		
		/// @name for BY_CELLS @{
		/// If some vertex has incident cells with different materials,
		/// the material with higher priority will be chosen.
		typedef int Priority; ///< must be >= 0
		typedef std::pair<Material, Priority> MaterialPriority;
		typedef size_t MaterialFlag;
		typedef std::map<MaterialFlag, MaterialPriority> MaterialMap;
		MaterialMap materialMap;
		/// @}
		
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
		///< legal number of outer characteristics
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

	struct Binary2DSeismograph {
		PhysicalQuantities::T quantityToWrite;
	} binary2DSeismograph;
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
	
	int dimensionality;
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
		std::vector<real> lengths;  ///< lengthes of cube in each direction
		std::vector<real> h;        ///< spatial steps in each direction
		std::vector<int>  sizes;    ///< number of nodes along each direction
		std::vector<real> startR;   ///< global coordinates of the first real node
		int borderSize = 0;         ///< number of virtual border nodes for one border point
		bool forceSequence = false; ///< behave such as Mpi::Size() == 1
	} cubicGrid;

	struct Cgal2DGrid {
		bool movable = false;       ///< deformable(true) or immutable(false) grid
		real spatialStep = 0;       ///< effective spatial step
		
		struct Body {
			typedef std::array<real, 2> Point;
			typedef std::vector<Point> Border;
			Body(const Border& outer_, const std::vector<Border>& inner_) :
					outer(outer_), inner(inner_) { }
			
			Border outer;               ///< outer border of the body
			std::vector<Border> inner;  ///< borders of the inner cavities
		};
		std::vector<Body> bodies; ///< list of bodies contained in the grid
	} cgal2DGrid;
	
	struct Cgal3DGrid {
		bool movable = false; ///< deformable(true) or immutable(false) grid
		
		/// @name Mesher properties @{
		enum class Mesher {
			CGAL_MESHER,
			INM_MESHER,
		} mesher;

		real spatialStep = 0; ///< effective spatial step FIXME InmMesher calculate spatial step!!!
		bool detectSharpEdges = false; ///< use true for figures with sharp edges
				///< (for CGAL mesher only)
		std::string fileName; ///< file with mesh to load from 
				///< (or with some initial data for perform meshing)
		real scale = 1; ///< denominator for given coordinates to scale the points
		/// @}
	} cgal3DGrid;

	/// list of statements to calculate on the same geometry
	std::vector<Statement> statements;
};

/** @} */
}

#endif // LIBGCM_TASK_HPP
