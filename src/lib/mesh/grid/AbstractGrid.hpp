#ifndef LIBGCM_ABSTRACTGRID_HPP
#define LIBGCM_ABSTRACTGRID_HPP


#include <lib/util/Logging.hpp>
#include <lib/util/task/Task.hpp>
#include <lib/Engine.hpp>

namespace gcm {

/**
 * Base class for all grids
 */
class AbstractGrid {
public:
	AbstractGrid(const Task&) { }
	virtual ~AbstractGrid() { }
};


}

#endif // LIBGCM_ABSTRACTGRID_HPP
