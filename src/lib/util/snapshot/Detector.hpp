#ifndef LIBGCM_DETECTOR_HPP
#define LIBGCM_DETECTOR_HPP

#include <lib/util/snapshot/Snapshotter.hpp>

namespace gcm {
	/**
	 * For inverse problem (aka zero-offset sensor, etc)
	 */
	template<class TMesh>
	class Detector : public Snapshotter {
	public:
		typedef typename TMesh::PdeVector       PdeVector;
		typedef typename TMesh::PartIterator    BorderIter;
		typedef typename TMesh::Model           Model;
		typedef typename Model::PdeVariables    PdeVariables;
		
		
		static const int DIMENSIONALITY = TMesh::DIMENSIONALITY;
		const std::string FILE_EXTENSION = std::string("bin");
		const std::string FOLDER_NAME    = std::string("1dseismo");

	protected:
		virtual void initializeImpl(const Task& task) override {
			assert_eq(task.detector.quantities.size(), 1); // more than one unsupported
			quantityToWrite = task.detector.quantities[0];
			detectionArea = task.detector.area;
		};
		virtual void beforeCalculationImpl(const Solver* solver) override {
			const TMesh* mesh = dynamic_cast<const TMesh*>(solver->getGrid());
			direction = DIMENSIONALITY - 1;
			indexOfDetectingSide = mesh->getSizes()(direction) - 1;
		};
		virtual void snapshotImpl(const AbstractGrid* mesh_, const int) override {
			const TMesh* mesh = dynamic_cast<const TMesh*>(mesh_);
			for (auto it = mesh->slice(direction, indexOfDetectingSide);
			          it != it.end(); ++it) {
				auto coords = mesh->coords(it);
				if (detectionArea->contains(coords)) {					
					detect(mesh->pde(it), coords);
				}
			}
			assert_ge(valuesInArea.size(), 1);
			real valueToWrite = std::accumulate(valuesInArea.begin(), valuesInArea.end(), 0)
					/ (real) valuesInArea.size();
			valuesInArea.clear();
			seismo.push_back((precision)valueToWrite);
		};
		virtual void afterCalculationImpl() override {
			FileUtils::openBinaryFileStream(fileStream, makeFileNameForSnapshot
					(seismoNumber, FILE_EXTENSION, FOLDER_NAME));
			FileUtils::writeStdVectorToBinaryFileStream(fileStream, seismo);
			FileUtils::closeFileStream(fileStream);
			seismoNumber++;
		};

	private:
		int seismoNumber = 0;
		std::vector<precision> seismo;
		std::vector<real> valuesInArea;

		std::shared_ptr<Area> detectionArea;
		int direction = 0;
		int indexOfDetectingSide = 0;
		
		PhysicalQuantities::T quantityToWrite;
		
		std::ofstream fileStream;
		std::ofstream fileStreamTest;
		
		void detect(const PdeVector pde, const linal::Vector3) {
			real realValue = PdeVariables::QUANTITIES.at(quantityToWrite).Get(pde);
			valuesInArea.push_back(realValue);
		};
	};
}

#endif // LIBGCM_DETECTOR_HPP
