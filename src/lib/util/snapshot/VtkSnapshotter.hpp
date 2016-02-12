#ifndef LIBGCM_VTKSNAPSHOTTER_HPP
#define LIBGCM_VTKSNAPSHOTTER_HPP

#include <lib/util/snapshot/Snapshotter.hpp>

#include <vtkSmartPointer.h>
#include <vtkPointData.h>


namespace gcm {
	class VtkSnapshotter : public Snapshotter {
	protected:
		static const int VtkVecSize = 3;
		std::vector<PhysicalQuantities::T> quantitiesToWrite = { };

		virtual void initializeImpl(const Task& task) override {
			quantitiesToWrite = task.quantitiesToWrite;
		};

		virtual void snapshotImpl(const Grid* _grid, const int step) = 0;

	};
}

#endif // LIBGCM_VTKSNAPSHOTTER_HPP
