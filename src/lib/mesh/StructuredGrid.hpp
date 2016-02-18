#ifndef LIBGCM_STRUCTUREDGRID_HPP
#define LIBGCM_STRUCTUREDGRID_HPP

#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>

#include <lib/mesh/AbstractGrid.hpp>

namespace gcm {
	/**
	 * Non-movable structured grid
	 */
	class StructuredGrid : public AbstractGrid {
	public:
		typedef vtkStructuredGrid             VtkGridType;
		typedef vtkXMLStructuredGridWriter    VtkWriterType;

		virtual ~StructuredGrid() { };
	protected:
		virtual void initializeImpl(const Task &task) = 0;
		virtual void beforeStageImpl() = 0;
		virtual void afterStageImpl() = 0;
		virtual void beforeStepImpl() = 0;
		virtual void afterStepImpl() = 0;

		virtual void recalculateMinimalSpatialStep() override { };

		virtual void recalculateMaximalLambda() = 0;
		virtual void applyInitialConditions(const Task& task) = 0;

		virtual void initializeImplImpl(const Task& task) = 0;

		USE_AND_INIT_LOGGER("gcm.StructuredGrid");
	};
}

#endif // LIBGCM_STRUCTUREDGRID_HPP
