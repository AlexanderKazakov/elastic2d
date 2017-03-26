#ifndef LIBGCM_ABSTRACTGRID_HPP
#define LIBGCM_ABSTRACTGRID_HPP


#include <libgcm/util/Enum.hpp>

namespace gcm {

/**
 * Base class for all grids
 */
class AbstractGrid {
public:
	
	/// unique number of the grid among others
	const GridId id;
	
	AbstractGrid(const GridId id_) : id(id_) { }
	virtual ~AbstractGrid() { }
	
};


}

#endif // LIBGCM_ABSTRACTGRID_HPP
