#ifndef LIBGCM_SLICESNAPSHOTTER_HPP
#define LIBGCM_SLICESNAPSHOTTER_HPP

#include <lib/util/snapshot/Snapshotter.hpp>

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

	const int DIMENSIONALITY = TMesh::DIMENSIONALITY;
	const std::string FILE_EXTENSION = std::string("txt");

protected:

	virtual void beforeStatementImpl(const Statement& statement) override {
		assert_eq(statement.detector.quantities.size(), 1); // more than one still unsupported
		quantityToWrite = statement.detector.quantities[0];
		detectionArea = statement.detector.area;
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
				std::to_string(mesh->id()),
				step, FILE_EXTENSION, "detector"), 
				std::vector<Values>({seismo}));
		
		// along Z-axis
		Int3 min = {mesh->sizes(0) / 2, mesh->sizes(1) / 2, 0};
		Int3 max = {1 + mesh->sizes(0) / 2, 1 + mesh->sizes(1) / 2, mesh->sizes(2)};
		std::vector<real> Vz, Szz, coordZ;
		for (auto it = mesh->box(min, max); it != it.end(); ++it) {
			coordZ.push_back(mesh->coords(it)(2));
			Vz.push_back(mesh->pdeVars(it).velocity(2));
			Szz.push_back(mesh->pdeVars(it).sigma(2, 2));
		}
		FileUtils::writeStdVectorsToTextFile(makeFileNameForSnapshot(
						std::to_string(mesh->id()),
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
