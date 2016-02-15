#ifndef LIBGCM_VTKSTRUCTUREDSNAPSHOTTER_HPP
#define LIBGCM_VTKSTRUCTUREDSNAPSHOTTER_HPP

#include <vtkStructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkXMLStructuredGridWriter.h>

#include <lib/util/snapshot/VtkSnapshotter.hpp>
#include <lib/util/Logging.hpp>

namespace gcm {
	template<class TGrid>
	class VtkStructuredSnapshotter : public VtkSnapshotter {
	protected:
		virtual void snapshotImpl(const Grid* _grid, const int step) override;

	private:
		USE_AND_INIT_LOGGER("gcm.VtkStructuredSnapshotter");
		const TGrid* grid;
		vtkSmartPointer<vtkStructuredGrid> vtkStrGrid;

		typedef void (VtkStructuredSnapshotter::*InsertFunc)(PhysicalQuantities::T quantity,
				vtkSmartPointer<vtkFloatArray> vtkArr, const typename TGrid::VtkIterator& it);

		void writeQuantity(PhysicalQuantities::T quantity, InsertFunc insertFunc, const int numOfComponents) {
			auto vtkArr = vtkSmartPointer<vtkFloatArray>::New();
			vtkArr->SetNumberOfComponents(numOfComponents);
			vtkArr->Allocate(linal::directProduct(grid->sizes), 0);
			vtkArr->SetName(PhysicalQuantities::NAME.at(quantity).c_str());
			for (auto it = grid->vtkBegin(); it != grid->vtkEnd(); ++it) {
				(this->*insertFunc)(quantity, vtkArr, it);
			}
			vtkStrGrid->GetPointData()->AddArray(vtkArr);
		};

		void insertVector(PhysicalQuantities::T quantity, vtkSmartPointer<vtkFloatArray> vtkArr,
		                  const typename TGrid::VtkIterator& it) {
			auto Get = TGrid::PdeVector::VECTORS.at(quantity).Get;
			auto linalVec = Get(grid->pde(it));
			float vtkVec[VtkVecSize];
			for (int i = 0; i < VtkVecSize; i++) {
				vtkVec[i] = (float) linalVec(i);
			}
			vtkArr->InsertNextTuple(vtkVec);
		};

		void insertQuantity(PhysicalQuantities::T quantity, vtkSmartPointer<vtkFloatArray> vtkArr,
		                    const typename TGrid::VtkIterator& it) {
			auto Get = TGrid::PdeVector::QUANTITIES.at(quantity).Get;
			vtkArr->InsertNextValue((float)Get(grid->pde(it)));
		};

		void insertOdeQuantity(PhysicalQuantities::T quantity, vtkSmartPointer<vtkFloatArray> vtkArr,
		                       const typename TGrid::VtkIterator& it) {
			auto Get = TGrid::Model::InternalOde::QUANTITIES.at(quantity).Get;
			vtkArr->InsertNextValue((float)Get(grid->ode(it)));
		};

		std::string makeFileNameForSnapshot(const int step);
	};
}

#endif // LIBGCM_VTKSTRUCTUREDSNAPSHOTTER_HPP
