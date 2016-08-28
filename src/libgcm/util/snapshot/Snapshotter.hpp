#ifndef LIBGCM_SNAPSHOTTER_HPP
#define LIBGCM_SNAPSHOTTER_HPP

#include <libgcm/util/infrastructure/infrastructure.hpp>
#include <libgcm/util/StringUtils.hpp>
#include <libgcm/util/FileUtils.hpp>
#include <libgcm/grid/AbstractGrid.hpp>
#include <libgcm/rheology/variables/GetSetter.hpp>
#include <libgcm/engine/GlobalVariables.hpp>
#include <libgcm/util/task/Task.hpp>


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
	const std::string SNAP = std::string("snap");
	const std::string SLASH = std::string("/");
	const std::string DOT = std::string(".");
	
	const int NUMBER_OF_DIGITS_IN_SNAP = 4;
	const int NUMBER_OF_DIGITS_IN_CORE = 2;
	///@}
	
	Snapshotter(const Task& task) {
		stepsPerSnap = task.globalSettings.stepsPerSnap;
	}
	virtual ~Snapshotter() { }
	
	void snapshot(const AbstractGrid* grid, const int step) {
		if (step % stepsPerSnap == 0) {
			snapshotImpl(grid, step);
		}
	}
	
	
protected:
	virtual void snapshotImpl(const AbstractGrid* grid, const int step) = 0;
	
	/** @return relative name for snapshot file */
	std::string makeFileNameForSnapshot(const std::string meshId, const int step,
			const std::string fileExtension, const std::string folder) const {
		
		std::string mesh = MESH + meshId;
		std::string snap = (step >= 0) ?
				SNAP + StringUtils::toString(step, NUMBER_OF_DIGITS_IN_SNAP) : "";
		std::string core = CORE + 
				StringUtils::toString(Mpi::Rank(), NUMBER_OF_DIGITS_IN_CORE);
		
		return SNAPSHOTS + SLASH + folder + SLASH +
				mesh + core + snap + DOT + fileExtension;
	}
	
	
private:
	int stepsPerSnap = 1;
	
};


}

#endif // LIBGCM_SNAPSHOTTER_HPP
