#ifndef LIBGCM_UNSTRUCTUREDGRID_HPP
#define LIBGCM_UNSTRUCTUREDGRID_HPP

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <lib/mesh/grid/AbstractGrid.hpp>
#include <lib/mesh/Elements.hpp>

namespace gcm {
/**
 * Movable unstructured grid
 */
class UnstructuredGrid : public AbstractGrid {
public:
	typedef vtkUnstructuredGrid          VtkGridType;
	typedef vtkXMLUnstructuredGridWriter VtkWriterType;

	/** 
	 * Simple size_t index iterator over all mesh nodes
	 */
	struct Iterator {
		size_t iter = 0;
		
		Iterator() : iter(0) { }
		Iterator(size_t value) : iter(value) { }
		
		const Iterator& operator*() const { return *this; }
		
		bool operator!=(const Iterator& other) const {
			return iter != other.iter;
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
	} Flags;

	UnstructuredGrid(const Task& task) : AbstractGrid(task) { }
	virtual ~UnstructuredGrid() { }

protected:
	USE_AND_INIT_LOGGER("gcm.UnstructuredGrid")
};


}

#endif // LIBGCM_UNSTRUCTUREDGRID_HPP
