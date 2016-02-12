#ifndef LIBGCM_BINARY2DSEISMOGRAPH_HPP
#define LIBGCM_BINARY2DSEISMOGRAPH_HPP

#include <fstream>

#include <lib/grid/StructuredGrid.hpp>
#include <lib/util/Logging.hpp>
#include <lib/util/snapshot/Snapshotter.hpp>

namespace gcm {
	/**
	 * For 2D binary seismography for inverse problem
	 */
	template<class TGrid>
	class Binary2DSeismograph : public Snapshotter {
	public:
		void startSeismo(const Task& task);
		void finishSeismo();

	protected:
		virtual void initializeImpl(const Task &task) override;
		virtual void snapshotImpl(const Grid* grid, const int step);

	private:
		USE_AND_INIT_LOGGER("gcm.Binary2DSeismograph");
		typedef float output_precision;

		std::ofstream seismoFileStream;

		output_precision* surface = nullptr; // values storage on current time step to write
		int sizeY = 0;
		output_precision hY = 0;
		output_precision tau = 0;
		int seismoNumber = 0;

		void writeHeadOfTable();
		void openSeismoFileStream(const std::string& fileName);
		void writeSurfaceToBuffer();
		void closeSeismoFileStream();
		std::string makeFileNameForSeismo();
	};
}

#endif // LIBGCM_BINARY2DSEISMOGRAPH_HPP
