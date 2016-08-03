#ifndef LIBGCM_STRUCTUREDGRID_HPP
#define LIBGCM_STRUCTUREDGRID_HPP

#include <lib/mesh/grid/AbstractGrid.hpp>

namespace gcm {

/**
 * Structured grid base
 */
class StructuredGrid : public AbstractGrid {
public:
	StructuredGrid(const Task& task) : AbstractGrid(task) { }
	virtual ~StructuredGrid() { }
	
};


}

#endif // LIBGCM_STRUCTUREDGRID_HPP
