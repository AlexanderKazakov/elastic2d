#ifndef LIBGCM_SLICESNAPSHOTTER_HPP
#define LIBGCM_SLICESNAPSHOTTER_HPP

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
		assert_eq(Mpi::Size() % 2, 1);
	}
	
	
	virtual void snapshotImpl(const AbstractGrid* mesh_, const int step) override {
		if (Mpi::Rank() != Mpi::Size() / 2) { return; }
		const TMesh* mesh = dynamic_cast<const TMesh*>(mesh_);
		assert_true(mesh);
		
		// upper detector
		int direction = DIMENSIONALITY - 1;
		int indexOfDetectingSide = mesh->sizes(direction) - 1;
		
		for (auto it = mesh->slice(direction, indexOfDetectingSide);
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
		seismo.push_back((precision)valueToWrite);
		
		FileUtils::writeStdVectorsToTextFile(makeFileNameForSnapshot(
				std::to_string(mesh->id),
				step, FILE_EXTENSION, "detector"), 
				std::vector<Values>({seismo}));
		
		// along Z-axis
		IntD min = mesh->sizes / 2;
		min(0) = 0;
		IntD max = mesh->sizes / 2 + IntD::Ones();
		max(DIMENSIONALITY) = mesh->sizes(DIMENSIONALITY);
		
		std::vector<real> Vz, Szz, coordZ;
		for (auto it = mesh->box(min, max); it != it.end(); ++it) {
			coordZ.push_back(mesh->coords(it)(DIMENSIONALITY));
			Vz.push_back(mesh->pdeVars(it).velocity(DIMENSIONALITY));
//			Szz.push_back(mesh->pdeVars(it).sigma(DIMENSIONALITY, DIMENSIONALITY));
		}
		FileUtils::writeStdVectorsToTextFile(makeFileNameForSnapshot(
						std::to_string(mesh->id),
						step, FILE_EXTENSION, "zaxis"), 
				std::vector<Values>({coordZ, Vz, Szz}));
	}
	
	
private:
	typedef std::vector<real> Values;
	Values seismo;
	Values valuesInArea;
	
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
