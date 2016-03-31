#ifndef LIBGCM_GETSETTER_HPP
#define LIBGCM_GETSETTER_HPP

#include <lib/linal/linal.hpp>

namespace gcm {

/**
 * @file
 * The structs below are made for classes like snapshotters, initial condition setters, etc.
 *
 * It provides unified for *all types* of Variables interface to information about their physical
 * quantities.
 * And so let one to get/set some quantity in an object of Variables without knowledge of its
 * structure:
 *
 * in every Variables class there is a static map { physical quantity <-> its GetSetter }
 *
 * The function pointers are supposed to reduce call overhead.
 * The map shouldn't be used at every access to every node, but just once before handling
 * a large portion of nodes.
 */

template<typename Variables>
struct GetSetter {
	typedef real (*Getter)(const Variables& variablesToGetFrom);
	typedef void (*Setter)(const real& value, Variables& variablesToSetTo);

	GetSetter(Getter _Get, Setter _Set) : Get(_Get), Set(_Set) { }

	Getter Get;
	Setter Set;
};


template<typename Variables>
struct Vector3GetSetter {
	typedef Real3 (*Getter)(const Variables& variablesToGetFrom);
	typedef void (*Setter)(const Real3& value, Variables& variablesToSetTo);

	Vector3GetSetter(Getter _Get, Setter _Set) : Get(_Get), Set(_Set) { }

	Getter Get;
	Setter Set;
};
}

#endif // LIBGCM_GETSETTER_HPP
