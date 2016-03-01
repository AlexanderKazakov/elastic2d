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
		
		
		const int DIMENSIONALITY = TMesh::DIMENSIONALITY;
		const std::string FILE_EXTENSION = std::string("bin");
		const std::string FOLDER_NAME    = std::string("1dseismo");

	protected:
		virtual void initializeImpl(const Task& task) override {
			direction = DIMENSIONALITY - 1;
			indexOfDetectingSide = task.sizes(direction) - 1;			
		};
		virtual void beforeStatementImpl(const Statement& statement) override {
			assert_eq(statement.detector.quantities.size(), 1); // more than one still unsupported
			quantityToWrite = statement.detector.quantities[0];
			detectionArea = statement.detector.area;
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
		virtual void afterStatementImpl() override {
			FileUtils::openBinaryFileStream(fileStream, makeFileNameForSnapshot
					(-1, FILE_EXTENSION, FOLDER_NAME));
			FileUtils::writeStdVectorToBinaryFileStream(fileStream, seismo);
			FileUtils::closeFileStream(fileStream);
		};

	private:
		std::vector<precision> seismo;
		std::vector<real> valuesInArea;

		std::shared_ptr<Area> detectionArea;
		int direction = 0;
		int indexOfDetectingSide = 0;
		
		PhysicalQuantities::T quantityToWrite;
		
		std::ofstream fileStream;
		
		void detect(const PdeVector pde, const linal::Vector3) {
			real realValue = PdeVariables::QUANTITIES.at(quantityToWrite).Get(pde);
			valuesInArea.push_back(realValue);
		};
	};
}

#endif // LIBGCM_DETECTOR_HPP
