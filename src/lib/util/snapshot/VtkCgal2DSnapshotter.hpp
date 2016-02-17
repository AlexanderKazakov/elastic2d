#ifndef LIBGCM_VTKCGAL2DSNAPSHOTTER_HPP
#define LIBGCM_VTKCGAL2DSNAPSHOTTER_HPP

#include <vtkXMLUnstructuredGridWriter.h>

#include <lib/util/snapshot/VtkSnapshotter.hpp>
#include <lib/util/Logging.hpp>

namespace gcm {
	template<class TGrid>
	class VtkCgal2DSnapshotter : public VtkSnapshotter {
	protected:
		virtual void snapshotImpl(const AbstractGrid* _grid, const int step) override;

	private:
		USE_AND_INIT_LOGGER("gcm.VtkCgal2DSnapshotter");
		const TGrid* grid;
		vtkSmartPointer<vtkUnstructuredGrid> vtkGrid;

		std::string makeFileNameForSnapshot(const int step);
	};
}

#endif // LIBGCM_VTKCGAL2DSNAPSHOTTER_HPP
