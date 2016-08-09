#ifndef LIBGCM_SNAPSHOTTER_HPP
#define LIBGCM_SNAPSHOTTER_HPP

#include <lib/util/StringUtils.hpp>
#include <lib/util/FileUtils.hpp>
#include <lib/mesh/grid/AbstractGrid.hpp>
#include <lib/rheology/variables/GetSetter.hpp>
#include <lib/GlobalVariables.hpp>

namespace gcm {

/**
 * Base class for snapshotters
 */

class Snapshotter {
public:
	typedef float precision;
	
	/** @name Details of snapshot file name */
	///@{
	const std::string MESH = std::string("mesh");
	const std::string CORE = std::string("core");
	const std::string SNAPSHOTS = std::string("snapshots");
	const std::string STATEMENT = std::string("statement");
	const std::string SNAP = std::string("snap");
	const std::string SLASH = std::string("/");
	const std::string DOT = std::string(".");
	
	const int NUMBER_OF_DIGITS_IN_SNAP = 4;
	const int NUMBER_OF_DIGITS_IN_CORE = 2;
	///@}
	
	
	void initialize(const Task& task) {
		if (task.statements.size() > 1) {
			writeStatementToFileName = true;
		}
	}
	
	
	void beforeStatement(const Statement& statement) {
		statementId = statement.id;
		stepsPerSnap = statement.globalSettings.stepsPerSnap;
		beforeStatementImpl(statement);
	}
	
	
	void snapshot(const AbstractGrid* grid, const int step) {
		if (step % stepsPerSnap == 0) {
			snapshotImpl(grid, step);
		}
	}
	
	
	virtual void afterStatement() { }
	virtual ~Snapshotter() { }
	
	
protected:
	virtual void beforeStatementImpl(const Statement&) { }
	virtual void snapshotImpl(const AbstractGrid* grid, const int step) = 0;
	
	/** @return relative name for snapshot file */
	std::string makeFileNameForSnapshot(const std::string meshId, const int step,
			const std::string fileExtension, const std::string folder) const {
		
		std::string mesh = MESH + meshId;
		std::string snap = (step >= 0) ?
				SNAP + StringUtils::toString(step, NUMBER_OF_DIGITS_IN_SNAP) : "";
		std::string core = CORE + 
				StringUtils::toString(Mpi::Rank(), NUMBER_OF_DIGITS_IN_CORE);
		std::string statement = writeStatementToFileName ?
				STATEMENT + statementId : "";
		
		return SNAPSHOTS + SLASH + folder + SLASH +
				mesh + core + statement + snap + DOT + fileExtension;
	}
	
	
private:
	bool writeStatementToFileName = false;
	std::string statementId = "";
	int stepsPerSnap = 1;
	
};


}

#endif // LIBGCM_SNAPSHOTTER_HPP
