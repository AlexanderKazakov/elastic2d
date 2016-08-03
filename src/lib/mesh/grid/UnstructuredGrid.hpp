#ifndef LIBGCM_UNSTRUCTUREDGRID_HPP
#define LIBGCM_UNSTRUCTUREDGRID_HPP

#include <lib/mesh/grid/AbstractGrid.hpp>
#include <lib/mesh/Elements.hpp>

namespace gcm {

/**
 * Unstructured grid base
 */
class UnstructuredGrid : public AbstractGrid {
public:
	
	/** 
	 * Simple plain size_t index iterator over all mesh nodes
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
	
	
	typedef union {
		unsigned char c = 0;
		struct {
			bool border;
		};
	} Flags; ///< currently unused
	
	
	UnstructuredGrid(const Task& task) : AbstractGrid(task) { }
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
