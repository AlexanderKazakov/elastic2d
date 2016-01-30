#ifndef LIBGCM_SNAPSHOTTER_HPP
#define LIBGCM_SNAPSHOTTER_HPP

#include <fstream>

#include <lib/grid/Grid.hpp>

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
		void snapshot(const Grid* grid, const int step) {
			if (enableSnapshotting && step % stepsPerSnap == 0) {
				snapshotImpl(grid, step);
			}
		};

		/**
		 * @param grid pointer to the grid for dump
		 * @param enableSnapshotting dump or do not dump snaps
		 */
		void initialize(const Task& task) {
			enableSnapshotting = task.enableSnapshotting;
			stepsPerSnap = task.stepsPerSnap;
		};
		virtual ~Snapshotter() { };

	protected:
		virtual void snapshotImpl(const Grid* grid, const int step) = 0;

	private:
		bool enableSnapshotting = false;
		int stepsPerSnap = 1;
	};
}

#endif // LIBGCM_SNAPSHOTTER_HPP
