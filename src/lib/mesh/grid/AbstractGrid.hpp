#ifndef LIBGCM_ABSTRACTGRID_HPP
#define LIBGCM_ABSTRACTGRID_HPP


#include <lib/util/Logging.hpp>
#include <lib/util/task/Task.hpp>


namespace gcm {

/**
 * Base class for all grids
 */
class AbstractGrid {
public:
	
	typedef size_t GridId;
	
	AbstractGrid(const Task& task,
					const Grids::T gridType_, const int dimensionality_,
					const GridId id_) :
			gridType(gridType_), dimensionality(dimensionality_), id(id_),
			modelId(determineModelType(task, id)),
			materialId(determineMaterialType(task, id)) { }
	
	virtual ~AbstractGrid() { }
	
	
	/// with this data one can cast mesh to its actual type
	const Grids::T gridType;
	const int dimensionality;
	
	/// unique number of the grid among others
	const GridId id;
	
	/// TODO - replace
	const Models::T modelId;
	const Materials::T materialId;
	
	
private:
	
	/// TODO - replace these duct tapes
	
	static Models::T determineModelType(const Task& task, const GridId id_) {
		Models::T ans = Models::T::ACOUSTIC;
		auto iter = task.bodies.find(id_);
		if (iter != task.bodies.end()) {
			ans = iter->second.modelId;
		}
		return ans;
	}
	static Materials::T determineMaterialType(const Task& task, const GridId id_) {
		Materials::T ans = Materials::T::ISOTROPIC;
		auto iter = task.bodies.find(id_);
		if (iter != task.bodies.end()) {
			ans = iter->second.materialId;
		}
		return ans;
	}
	
	
};


}


#endif // LIBGCM_ABSTRACTGRID_HPP
