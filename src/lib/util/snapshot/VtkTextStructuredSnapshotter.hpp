#ifndef LIBGCM_VTKTEXTSTRUCTUREDSNAPSHOTTER_HPP
#define LIBGCM_VTKTEXTSTRUCTUREDSNAPSHOTTER_HPP

#include <fstream>

#include <lib/grid/StructuredGrid.hpp>
#include <lib/util/snapshot/Snapshotter.hpp>
#include <lib/util/Logging.hpp>

namespace gcm {
	template<class TNode> class StructuredGrid;

	template<class TNode>
	class VtkTextStructuredSnapshotter : public Snapshotter<TNode> {
	public:
		/**
		 * Write vtk snapshot. The simplest but not the most efficient
		 * snapshotting performed by writing data to file as a text.
		 * VTK library is not necessary for such snapshotting.
		 */
		virtual void snapshotImpl(const int step) override;

	private:
		USE_AND_INIT_LOGGER("gcm.VtkTextStructuredSnapshotter");
		StructuredGrid<TNode>* sGrid;

		void writeQuantity(const std::string name, const typename GetSetter<typename TNode::Variables>::Getter Get);

		std::fstream snapshotFileStream;
		void openSnapshotFileStream(const std::string& fileName);
		void closeSnapshotFileStream();

		std::string makeFileNameForSnapshot(const int step);
	};
}

#endif // LIBGCM_VTKTEXTSTRUCTUREDSNAPSHOTTER_HPP
