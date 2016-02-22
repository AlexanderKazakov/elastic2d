#ifndef LIBGCM_SNAPSHOTTER_HPP
#define LIBGCM_SNAPSHOTTER_HPP

#include <fstream>

#include <lib/mesh/grid/AbstractGrid.hpp>

namespace gcm {
	/**
	 * Base class for snapshotters
	 */
	class Snapshotter {
	public:
		/**
		 * Write snapshot for specified time step
		 * @param step number of time step
		 */
		void snapshot(const AbstractGrid* grid, const int step) {
			if (enableSnapshotting && step % stepsPerSnap == 0) {
				snapshotImpl(grid, step);
			}
		};

		void initialize(const Task& task) {
			enableSnapshotting = task.enableSnapshotting;
			stepsPerSnap = task.stepsPerSnap;
			initializeImpl(task);
		};
		virtual ~Snapshotter() { };

	protected:
		virtual void snapshotImpl(const AbstractGrid* grid, const int step) = 0;
		virtual void initializeImpl(const Task& task) = 0;

	private:
		bool enableSnapshotting = false;
		int stepsPerSnap = 1;
	};
}

#endif // LIBGCM_SNAPSHOTTER_HPP
