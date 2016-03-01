#ifndef LIBGCM_SNAPSHOTTER_HPP
#define LIBGCM_SNAPSHOTTER_HPP

#include <lib/util/StringUtils.hpp>
#include <lib/util/FileUtils.hpp>
#include <lib/mesh/grid/AbstractGrid.hpp>
#include <lib/numeric/solvers/Solver.hpp>

namespace gcm {
	/**
	 * Base class for snapshotters
	 */
	class Snapshotter {
	public:
		typedef float precision;
		// details of snapshot file name
		const std::string CORE      = std::string("core");
		const std::string SNAPSHOTS = std::string("snapshots");
		const std::string STEP      = std::string("step");
		const std::string SLASH     = std::string("/");
		const std::string DOT       = std::string(".");
		
		void initialize(const Task& task) {
			enableSnapshotting = task.enableSnapshotting;
			stepsPerSnap = task.stepsPerSnap;
			if (enableSnapshotting) initializeImpl(task);
		};
		void beforeCalculation(const Solver* solver) {
			if (enableSnapshotting) beforeCalculationImpl(solver);
		};
		void snapshot(const AbstractGrid* grid, const int step) {
			if (enableSnapshotting && step % stepsPerSnap == 0) {
				snapshotImpl(grid, step);
			}
		};
		void afterCalculation() {
			if (enableSnapshotting) afterCalculationImpl();
		};
		virtual ~Snapshotter() { };

	protected:
		virtual void initializeImpl(const Task&) { };
		virtual void beforeCalculationImpl(const Solver*) { };
		virtual void snapshotImpl(const AbstractGrid* grid, const int step) = 0;
		virtual void afterCalculationImpl() { };

		/** @return relative name for snapshot file */
		std::string makeFileNameForSnapshot(const int step_, 
				const std::string fileExtension, const std::string folder) const {
			const auto step = StringUtils::toString(step_, 4);
			const auto core = StringUtils::toString(MPI::COMM_WORLD.Get_rank(), 2);
			return SNAPSHOTS + SLASH + folder + SLASH + CORE + core +
			       STEP + step + DOT + fileExtension;
		};
		
	private:
		bool enableSnapshotting = false;
		int stepsPerSnap = 1;
	};
}

#endif // LIBGCM_SNAPSHOTTER_HPP
