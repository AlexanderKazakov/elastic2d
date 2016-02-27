#ifndef LIBGCM_CONCEPTS_HPP
#define LIBGCM_CONCEPTS_HPP

#include <map>
#include <string>

/**
 * The file contains structs with enumerations and some properties for physical concepts used by the program.
 * The aim is to formalize and unify interconnection between parts of the program.
 * Use it every time when parts of the program communicate to each other by this concepts.
 * Enum class called T (type), map with string names called NAME.
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

			SIZE_OF_ENUM /* This MUST be at the last position here! DO NOT insert anything after that!!! */
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

			SIZE_OF_ENUM /* This MUST be at the last position here! DO NOT insert anything after that!!! */
		};

		/** string names of concepts */
		static const std::map<T, std::string> NAME;
	};
	
	enum class DIRECTION {X = 0, Y = 1, Z = 2};
};


#endif // LIBGCM_CONCEPTS_HPP
