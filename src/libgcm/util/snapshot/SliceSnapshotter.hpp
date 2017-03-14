#ifndef LIBGCM_SLICESNAPSHOTTER_HPP
#define LIBGCM_SLICESNAPSHOTTER_HPP

#include <numeric>

#include <libgcm/util/snapshot/Snapshotter.hpp>

namespace gcm {

/**
 * Write Szz and Vz along Z-axis into one text file
 * and mean detector value from the top of cube (along Z) into another one.
 */
template<class TMesh>
class SliceSnapshotter : public Snapshotter {
public:
	typedef typename TMesh::PdeVector    PdeVector;
	typedef typename TMesh::Model        Model;
	typedef typename Model::PdeVariables PdeVariables;
	typedef typename TMesh::Grid::IntD   IntD;
	
	const int DIMENSIONALITY = TMesh::DIMENSIONALITY;
	const std::string FILE_EXTENSION = std::string("txt");
	
	
	SliceSnapshotter(const Task& task) : Snapshotter(task) {
		assert_eq(task.detector.quantities.size(), 1); // more than one still unsupported
		quantityToWrite = task.detector.quantities[0];
		detectionArea = task.detector.area;
		seismo.clear();
		gridId = task.detector.gridId;
		assert_eq(Mpi::Size() % 2, 1);
	}
	
	
	virtual void snapshotImpl(const AbstractGrid* mesh_, const int step) override {
		if (Mpi::Rank() != Mpi::Size() / 2) { return; }
		const TMesh* mesh = dynamic_cast<const TMesh*>(mesh_);
		assert_true(mesh);
		int direction = DIMENSIONALITY - 1;
		
		
		// along Z-axis
		IntD min = mesh->sizes / 2;
		min(direction) = 0;
		IntD max = min + IntD::Ones();
		max(direction) = mesh->sizes(direction);
		
		std::vector<real> Vz, /*Szz,*/ coordZ;
		for (auto it = mesh->box(min, max); it != it.end(); ++it) {
			coordZ.push_back(mesh->coords(it)(direction));
			Vz.push_back(mesh->pdeVars(it).velocity(direction));
		}
		FileUtils::writeStdVectorsToTextFile(makeFileNameForSnapshot(
						std::to_string(mesh->id),
						step, FILE_EXTENSION, "zaxis"), 
				std::vector<Values>({coordZ, Vz/*, Szz*/}));
		
		
		// upper detector
		if (mesh_->id != gridId) { return; }
		for (auto it = mesh->rightBorder(direction);
				it != it.end(); ++it) {
			auto coords = mesh->coords(it);
			if (detectionArea->contains(coords)) {
				detect(mesh->pde(it), coords);
			}
		}
		assert_ge(valuesInArea.size(), 1);
		real valueToWrite =
				std::accumulate(valuesInArea.begin(), valuesInArea.end(), 0.0) /
				(real) valuesInArea.size();
		valuesInArea.clear();
		times.push_back(Clock::Time());
		seismo.push_back((precision)valueToWrite);
		
		FileUtils::writeStdVectorsToTextFile(makeFileNameForSnapshot(
				std::to_string(mesh->id),
				step, FILE_EXTENSION, "detector"), 
				std::vector<Values>({times, seismo}));
		
	}
	
	
private:
	typedef std::vector<real> Values;
	size_t gridId;
	Values times; //< array of time moments
	Values seismo; //< array of values at time moments
	Values valuesInArea; //< helper
	
	std::shared_ptr<Area> detectionArea;
	
	PhysicalQuantities::T quantityToWrite;
	
	std::ofstream fileStream;
	
	void detect(const PdeVector pde, const Real3) {
		real realValue = PdeVariables::QUANTITIES.at(quantityToWrite).Get(pde);
		valuesInArea.push_back(realValue);
	}
	
};


}

#endif // LIBGCM_SLICESNAPSHOTTER_HPP
