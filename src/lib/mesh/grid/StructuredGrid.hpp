#ifndef LIBGCM_STRUCTUREDGRID_HPP
#define LIBGCM_STRUCTUREDGRID_HPP

#include <lib/mesh/grid/AbstractGrid.hpp>

namespace gcm {
/**
 * Non-movable structured grid
 */
class StructuredGrid : public AbstractGrid {
public:
	StructuredGrid(const Task& task) : AbstractGrid(task) { }
	virtual ~StructuredGrid() { }

protected:
	USE_AND_INIT_LOGGER("gcm.StructuredGrid")
};


}

#endif // LIBGCM_STRUCTUREDGRID_HPP
