#ifndef LIBGCM_CONCEPTS_HPP
#define LIBGCM_CONCEPTS_HPP


#include <map>
#include <string>


/**
 * @file
 * The file contains structs with enumerations and some properties for physical concepts used by the program.
 * The aim is to formalize and unify interconnection between parts of the program.
 * Use it every time when parts of the program communicate to each other by this concepts.
 * Enum class called T (type), map with string names called NAME.
 * TODO - improve
 */

namespace gcm {

/**
 * For all physical quantities used in the program.
 */
struct PhysicalQuantities {
	/** Type */
	enum class T {
		VELOCITY,
		FORCE,
		
		Vx /* Velocity component along x-axis */,
		Vy /* Velocity component along y-axis */,
		Vz /* Velocity component along z-axis */,
		
		Sxx /* Component of tension tensor */,
		Sxy /* Component of tension tensor */,
		Sxz /* Component of tension tensor */,
		Syy /* Component of tension tensor */,
		Syz /* Component of tension tensor */,
		Szz /* Component of tension tensor */,
		
		RHO /* Density */,
		PRESSURE /* Pressure = -Sxx (1D), -(Sxx + Syy)/2 (2D), -(Sxx + Syy + Szz)/3 (3D) */,
		
		DAMAGE_MEASURE /* from continual damage model */,
		
	};
	/** string names of concepts */
	static const std::map<T, std::string> NAME;
};

/**
 * For all types of waves used in the program.
 */
struct Waves {
	/** Type */
	enum class T {
		P_FORWARD, P_BACKWARD,
		S1_FORWARD, S1_BACKWARD,
		S2_FORWARD, S2_BACKWARD,
		
	};
	/** string names of concepts */
	static const std::map<T, std::string> NAME;
};


struct ContactConditions {
	/** Type */
	enum class T {
		ADHESION,
		SLIDE
	};
//	/** string names of concepts */
//	static const std::map<T, std::string> NAME;
};

struct BorderConditions {
	/** Type */
	enum class T {
		FIXED_FORCE,
		FIXED_VELOCITY
	};
	/** string names of concepts */
	static const std::map<T, std::string> NAME;
};

/**
 * For all types of materials used in the program.
 */
struct Materials {
	/** Type */
	enum class T {
		ISOTROPIC,
		ORTHOTROPIC,
		
	};
	/** string names of concepts */
	static const std::map<T, std::string> NAME;
};

/**
 * For all types of grids used in the program.
 */
struct Grids {
	/** Type */
	enum class T {
		CUBIC,
		SIMPLEX,
		
	};
	/** string names of concepts */
	static const std::map<T, std::string> NAME;
};

/**
 * For all types of rheology models used in the program.
 */
struct Models {
	/** Type */
	enum class T {
		ELASTIC,
		ACOUSTIC,
		MAXWELL_ACOUSTIC,
		
	};
	/** string names of concepts */
	static const std::map<T, std::string> NAME;
};

/**
 * For all types of snapshotters used in the program.
 */
struct Snapshotters {
	/** Type */
	enum class T {
		VTK,
		DETECTOR,
		SLICESNAP,
		
	};

	/** string names of concepts */
	static const std::map<T, std::string> NAME;
};


/**
 * For all types of ODE's and correctors used in the program.
 */
struct Odes {
	/** Type */
	enum class T {
		MAXWELL_VISCOSITY,
		CONTINUAL_DAMAGE,
		IDEAL_PLASTIC_FLOW,
		
	};
	
//	/** string names of concepts */
//	static const std::map<T, std::string> NAME;
};


/**
 * Different border-contact nodes calculation modes for simplex engine
 */
enum class BorderCalcMode {
	LOCAL_BASIS,
	GLOBAL_BASIS,
};


/**
 * Different types of gcm-method to use for calculations
 */
enum class GcmType {
	ADVECT_RIEMANN_INVARIANTS,
	ADVECT_PDE_VECTORS,
};


/**
 * Splitting by directions approach
 */
enum class SplittingType {
	SUMM,   ///< (A_1 * u + A_2 * u) / 2
	PRODUCT ///< A_2 * (A_1 * u)
};

}


#endif // LIBGCM_CONCEPTS_HPP
