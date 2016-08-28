#ifndef LIBGCM_UNSTRUCTUREDGRID_HPP
#define LIBGCM_UNSTRUCTUREDGRID_HPP

#include <libgcm/grid/AbstractGrid.hpp>
#include <libgcm/util/Elements.hpp>


namespace gcm {

/**
 * Unstructured grids base
 */
class UnstructuredGrid : public AbstractGrid {
public:
	
	/** 
	 * Simple plain size_t index iterator over all grid nodes
	 */
	struct Iterator {
		typedef size_t Index;
		
		Index iter = 0;
		
		Iterator(size_t value = 0) : iter(value) { }
		
		operator Index() const { return iter; }
		
		const Iterator& operator*() const { return *this; }
		
		bool operator==(const Iterator& other) const {
			return iter == other.iter;
		}
		
		bool operator!=(const Iterator& other) const {
			return !( (*this) == other );
		}
		
		bool operator<(const Iterator& other) const {
			return iter < other.iter;
		}
	
		Iterator& operator++() {
			iter++;
			return (*this);
		}
	};
	
	
	UnstructuredGrid(const GridId id_) : AbstractGrid(id_) { }
	virtual ~UnstructuredGrid() { }
	
	
	/**
	 * @param it begin() <= iterator < end()
	 * @return index of node by the Iterator in std::vector storage
	 */
	size_t getIndex(const Iterator& it) const {
		return it.iter;
	}
	
};


}

#endif // LIBGCM_UNSTRUCTUREDGRID_HPP
