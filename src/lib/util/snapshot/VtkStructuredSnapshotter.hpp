#ifndef LIBGCM_VTKSTRUCTUREDSNAPSHOTTER_HPP
#define LIBGCM_VTKSTRUCTUREDSNAPSHOTTER_HPP

#include <vtkStructuredGrid.h>

#include <lib/util/snapshot/Snapshotter.hpp>
#include <lib/util/Logging.hpp>
#include <vtkSmartPointer.h>

namespace gcm {
	template<class TGrid>
	class VtkStructuredSnapshotter : public Snapshotter {
	protected:
		virtual void snapshotImpl(const Grid* _grid, const int step) override;

	private:
		USE_AND_INIT_LOGGER("gcm.VtkStructuredSnapshotter");
		std::vector<PhysicalQuantities::T> quantitiesToWrite = {PhysicalQuantities::T::PRESSURE,
		PhysicalQuantities::T::Sxx, PhysicalQuantities::T::Sxy, PhysicalQuantities::T::Syy};
		const TGrid* grid;
		vtkSmartPointer<vtkStructuredGrid> vtkStrGrid;
		static const int vtkVecSize = 3;

		void writeVector(const std::string name,
		                 const typename Vector3GetSetter<typename TGrid::Model::Variables>::Getter Get);
		void writeQuantity(const std::string name,
		                   const typename GetSetter<typename TGrid::Model::Variables>::Getter Get);

		void writeQuantity(const std::string name,
		            const typename GetSetter<typename TGrid::OdeVariables>::Getter Get);

		std::string makeFileNameForSnapshot(const int step);
	};
}

#endif // LIBGCM_VTKSTRUCTUREDSNAPSHOTTER_HPP
