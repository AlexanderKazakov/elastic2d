#ifndef LIBGCM_VTKTEXTSTRUCTUREDSNAPSHOTTER_HPP
#define LIBGCM_VTKTEXTSTRUCTUREDSNAPSHOTTER_HPP

#include <fstream>
#include <mpi.h>

#include "lib/grid/StructuredGrid.hpp"
#include "lib/snapshot/Snapshotter.hpp"
#include "lib/util/Logging.hpp"

namespace gcm {
	template<class TModel> class StructuredGrid;

	template<class TModel>
	class VtkTextStructuredSnapshotter : public Snapshotter {
	public:
		/**
		 * Write vtk snapshot. The simplest but not the most efficient
		 * snapshotting performed by writing data to file as a text.
		 * VTK library is not needed for such snapshotting.
		 */
		virtual void snapshotImpl(const std::string& fileName) override;

	private:
		USE_AND_INIT_LOGGER("gcm.VtkTextStructuredSnapshotter");
		StructuredGrid<TModel>* structuredGrid;

		void writeScalar(const std::string name, const int index);

		static const int VTK_VECTOR_SIZE = 3;
		void writeVector(const std::string& name, const int index, const int size);

		std::fstream snapshotFileStream;
		void openSnapshotFileStream(const std::string& fileName);
		void closeSnapshotFileStream();
	};
}

#endif // LIBGCM_VTKTEXTSTRUCTUREDSNAPSHOTTER_HPP
