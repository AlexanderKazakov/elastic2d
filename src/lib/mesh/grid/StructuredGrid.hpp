#ifndef LIBGCM_STRUCTUREDGRID_HPP
#define LIBGCM_STRUCTUREDGRID_HPP

#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>

#include <lib/mesh/grid/AbstractGrid.hpp>

namespace gcm {
	/**
	 * Non-movable structured grid
	 */
	class StructuredGrid : public AbstractGrid {
	public:
		typedef vtkStructuredGrid             VtkGridType;
		typedef vtkXMLStructuredGridWriter    VtkWriterType;

		StructuredGrid(const Task& task) : AbstractGrid(task) { }
		virtual ~StructuredGrid() { }

	protected:
		virtual void recalculateMinimalSpatialStep() override { }

		USE_AND_INIT_LOGGER("gcm.StructuredGrid")
	};
}

#endif // LIBGCM_STRUCTUREDGRID_HPP
