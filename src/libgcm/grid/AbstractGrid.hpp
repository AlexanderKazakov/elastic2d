#ifndef LIBGCM_ABSTRACTGRID_HPP
#define LIBGCM_ABSTRACTGRID_HPP


namespace gcm {

/**
 * Base class for all grids
 */
class AbstractGrid {
public:
	
	/// unique number of the grid among others
	typedef size_t GridId;
	const GridId id;
	
	AbstractGrid(const GridId id_) : id(id_) { }
	virtual ~AbstractGrid() { }
	
};


}

#endif // LIBGCM_ABSTRACTGRID_HPP
