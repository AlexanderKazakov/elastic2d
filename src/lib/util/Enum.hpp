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
		BIN2DSEISM,
		DETECTOR,
		SLICESNAP,
		
	};

	/** string names of concepts */
	static const std::map<T, std::string> NAME;
};

/**
 * For all types of solvers used in the program.
 */
struct Solvers {
	/** Type */
	enum class T {
		DEFAULT,
	};
};

/**
 * For all types of meshes used in the program.
 */
struct Meshes {
	/** Type */
	enum class T {
		DEFAULT,
	};
};

enum class DIRECTION {X = 0, Y = 1, Z = 2};


}


#endif // LIBGCM_CONCEPTS_HPP
