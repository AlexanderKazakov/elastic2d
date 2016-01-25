#ifndef LIBGCM_SNAPSHOTTER_HPP
#define LIBGCM_SNAPSHOTTER_HPP

#include <fstream>

#include <lib/grid/Grid.hpp>

namespace gcm {
	/**
	 * Base class for snapshotters
	 */
	template <class TModel>
	class Snapshotter {
	public:
		/**
		 * Write snapshot for specified time step
		 * @param step number of time step
		 */
		void snapshot(const int step);

		/**
		 * @param grid pointer to the grid for dump
		 * @param enableSnapshotting dump or do not dump snaps
		 */
		void initialize(Grid<TModel>* _grid, bool _enableSnapshotting);
		virtual ~Snapshotter() { };

	protected:
		Grid<TModel>* grid = nullptr;
		virtual void snapshotImpl(const std::string& fileName) = 0;

	private:
		bool enableSnapshotting = false;
		std::string makeFileNameForSnapshot(const int step);
	};
}

#endif // LIBGCM_SNAPSHOTTER_HPP
