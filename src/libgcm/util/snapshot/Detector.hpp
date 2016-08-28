#ifndef LIBGCM_DETECTOR_HPP
#define LIBGCM_DETECTOR_HPP

#include <libgcm/util/snapshot/Snapshotter.hpp>

namespace gcm {

/**
 * For inverse problem (aka zero-offset sensor, etc)
 */
template<class TMesh>
class Detector : public Snapshotter {
public:
	typedef typename TMesh::PdeVector    PdeVector;
	typedef typename TMesh::Model        Model;
	typedef typename Model::PdeVariables PdeVariables;
	
	const int DIMENSIONALITY = TMesh::DIMENSIONALITY;
	const std::string FILE_EXTENSION = std::string("bin");
	const std::string FOLDER_NAME = std::string("1dseismo");
	
	
	Detector(const Task& task) : Snapshotter(task) {
		assert_eq(task.detector.quantities.size(), 1); // more than one still unsupported
		quantityToWrite = task.detector.quantities[0];
		detectionArea = task.detector.area;
	}
	
	virtual void snapshotImpl(const AbstractGrid* mesh_, const int step) override {
		const TMesh* mesh = dynamic_cast<const TMesh*>(mesh_);
		assert_true(mesh);
		
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
		real valueToWrite = std::accumulate(
				valuesInArea.begin(), valuesInArea.end(), 0.0) /
						(real) valuesInArea.size();
		valuesInArea.clear();
		seismo.push_back((precision)valueToWrite);
		
		FileUtils::openTextFileStream(fileStream, makeFileNameForSnapshot
				("?", step, FILE_EXTENSION, FOLDER_NAME));
		FileUtils::writeStdVectorToTextFileStream(fileStream, seismo);
		FileUtils::closeFileStream(fileStream);
	}
	
	
private:
	std::vector<precision> seismo;
	std::vector<real> valuesInArea;
	
	std::shared_ptr<Area> detectionArea;
	
	PhysicalQuantities::T quantityToWrite;
	
	std::ofstream fileStream;
	
	void detect(const PdeVector pde, const Real3) {
		real realValue = PdeVariables::QUANTITIES.at(quantityToWrite).Get(pde);
		valuesInArea.push_back(realValue);
	}
	
};


}

#endif // LIBGCM_DETECTOR_HPP
