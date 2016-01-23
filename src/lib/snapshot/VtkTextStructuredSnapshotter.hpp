#ifndef LIBGCM_VTKTEXTSTRUCTUREDSNAPSHOTTER_HPP
#define LIBGCM_VTKTEXTSTRUCTUREDSNAPSHOTTER_HPP

#include <fstream>

#include <lib/grid/StructuredGrid.hpp>
#include <lib/snapshot/Snapshotter.hpp>
#include <lib/util/Logging.hpp>

namespace gcm {
	template<class TModel> class StructuredGrid;

	template<class TModel>
	class VtkTextStructuredSnapshotter : public Snapshotter<TModel> {
	public:
		/**
		 * Write vtk snapshot. The simplest but not the most efficient
		 * snapshotting performed by writing data to file as a text.
		 * VTK library is not necessary for such snapshotting.
		 */
		virtual void snapshotImpl(const std::string& fileName) override;

	private:
		USE_AND_INIT_LOGGER("gcm.VtkTextStructuredSnapshotter");
		StructuredGrid<TModel>* sGrid;

		void writeQuantity(const std::string name, const typename TModel::Getter Get);

		std::fstream snapshotFileStream;
		void openSnapshotFileStream(const std::string& fileName);
		void closeSnapshotFileStream();
	};
}

#endif // LIBGCM_VTKTEXTSTRUCTUREDSNAPSHOTTER_HPP
