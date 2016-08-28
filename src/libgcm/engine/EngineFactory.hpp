#ifndef LIBGCM_ENGINEFACTORY_HPP
#define LIBGCM_ENGINEFACTORY_HPP

#include <libgcm/engine/cubic/Engine.hpp>
#include <libgcm/engine/simplex/Engine.hpp>
#include <libgcm/grid/simplex/cgal/CgalTriangulation.hpp>


namespace gcm {

inline std::shared_ptr<AbstractEngine> createEngine(const Task& task) {
	
	switch (task.globalSettings.gridId) {
	
	case Grids::T::CUBIC:
	switch (task.globalSettings.dimensionality) {
		case 1: return std::make_shared<cubic::Engine<1>>(task);
		case 2: return std::make_shared<cubic::Engine<2>>(task);
		case 3: return std::make_shared<cubic::Engine<3>>(task);
		default: THROW_INVALID_ARG("Invalid space dimensionality");
	}
	
	case Grids::T::SIMPLEX:
	switch (task.globalSettings.dimensionality) {
		case 1: THROW_UNSUPPORTED("Unsupported space dimensionality");
		case 2: return std::make_shared<simplex::Engine<2, CgalTriangulation>>(task);
		case 3: return std::make_shared<simplex::Engine<3, CgalTriangulation>>(task);
		default: THROW_INVALID_ARG("Invalid space dimensionality");
	}
	
	default:
		THROW_UNSUPPORTED("Unknown type of grid");
	
	}
	
}


}

#endif // LIBGCM_ENGINEFACTORY_HPP
