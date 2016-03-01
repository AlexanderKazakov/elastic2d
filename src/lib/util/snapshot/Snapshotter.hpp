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
		const std::string STATEMENT = std::string("statement");
		const std::string SNAP      = std::string("snap");
		const std::string SLASH     = std::string("/");
		const std::string DOT       = std::string(".");
		
		const int NUMBER_OF_DIGITS_IN_SNAP = 4;
		const int NUMBER_OF_DIGITS_IN_CORE = 2;
		
		void initialize(const Task& task) {
			enableSnapshotting = task.enableSnapshotting;
			if (task.statements.size() > 1) writeStatementToFileName = true;
			if (enableSnapshotting) initializeImpl(task);
		};
		void beforeStatement(const Statement& statement) {
			statementId = statement.id;
			stepsPerSnap = statement.stepsPerSnap;
			if (enableSnapshotting) beforeStatementImpl(statement);
		};
		void snapshot(const AbstractGrid* grid, const int step) {
			if (enableSnapshotting && step % stepsPerSnap == 0) {
				snapshotImpl(grid, step);
			}
		};
		void afterStatement() {
			if (enableSnapshotting) afterStatementImpl();
		};
		virtual ~Snapshotter() { };

	protected:
		virtual void initializeImpl(const Task&) { };
		virtual void beforeStatementImpl(const Statement&) { };
		virtual void snapshotImpl(const AbstractGrid* grid, const int step) = 0;
		virtual void afterStatementImpl() { };

		/** @return relative name for snapshot file */
		std::string makeFileNameForSnapshot(const int step, 
				const std::string fileExtension, const std::string folder) const {
			
			std::string snap = (step >= 0) ? 
				SNAP + StringUtils::toString(step, NUMBER_OF_DIGITS_IN_SNAP) : "";
			std::string core = CORE + 
				StringUtils::toString(MPI::COMM_WORLD.Get_rank(), NUMBER_OF_DIGITS_IN_CORE);
			std::string statement = writeStatementToFileName ?
				STATEMENT + statementId : "";
			
			return SNAPSHOTS + SLASH + folder + SLASH + core + statement + snap + DOT + fileExtension;
		};
		
	private:
		bool enableSnapshotting = false;
		bool writeStatementToFileName = false;
		std::string statementId = "";
		int stepsPerSnap = 1;
	};
}

#endif // LIBGCM_SNAPSHOTTER_HPP
