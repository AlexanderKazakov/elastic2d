#ifndef LIBGCM_UNSTRUCTUREDGRID_HPP
#define LIBGCM_UNSTRUCTUREDGRID_HPP

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <lib/mesh/AbstractGrid.hpp>

namespace gcm {
	/**
	 * Movable unstructured grid
	 */
	class UnstructuredGrid : public AbstractGrid {
	public:
		typedef vtkUnstructuredGrid             VtkGridType;
		typedef vtkXMLUnstructuredGridWriter    VtkWriterType;

		struct Iterator {
			size_t iter = 0;
			Iterator(size_t value) : iter(value) { };
			const Iterator& operator*() { return *this; };
			bool operator!=(const Iterator& other) const {
				return iter != other.iter;
			};
			Iterator& operator++() {
				iter++;
				return (*this);
			};
		};

		virtual ~UnstructuredGrid() { };

	protected:
		virtual void initializeImpl(const Task &task) = 0;
		virtual void beforeStageImpl() = 0;
		virtual void afterStageImpl() = 0;
		virtual void beforeStepImpl() = 0;
		virtual void afterStepImpl() = 0;

		virtual void recalculateMinimalSpatialStep() = 0;

		virtual void recalculateMaximalLambda() = 0;
		virtual void applyInitialConditions(const Task& task) = 0;

		virtual void initializeImplImpl(const Task& task) = 0;

		USE_AND_INIT_LOGGER("gcm.UnstructuredGrid");
	};
}

#endif // LIBGCM_UNSTRUCTUREDGRID_HPP
