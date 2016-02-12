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

		typedef void (VtkStructuredSnapshotter::*InsertFunc)(PhysicalQuantities::T quantity, vtkSmartPointer<vtkFloatArray> vtkArr,
		                                                     const int x, const int y, const int z);

		void writeQuantity(PhysicalQuantities::T quantity, InsertFunc insertFunc, const int numOfComponents) {
			auto vtkArr = vtkSmartPointer<vtkFloatArray>::New();
			vtkArr->SetNumberOfComponents(numOfComponents);
			vtkArr->Allocate(linal::directProduct(grid->sizes), 0);
			vtkArr->SetName(PhysicalQuantities::NAME.at(quantity).c_str());

			for (int z = 0; z < grid->sizes(2); z++) {
				for (int y = 0; y < grid->sizes(1); y++) {
					for (int x = 0; x < grid->sizes(0); x++) {
						(this->*insertFunc)(quantity, vtkArr, x, y, z);
					}
				}
			}
			vtkStrGrid->GetPointData()->AddArray(vtkArr);
		};

		void insertVector(PhysicalQuantities::T quantity, vtkSmartPointer<vtkFloatArray> vtkArr,
		                         const int x, const int y, const int z) {
			auto Get = TGrid::PdeVector::VECTORS.at(quantity).Get;
			auto linalVec = Get(grid->get(x, y, z));
			float vtkVec[VtkVecSize];
			for (int i = 0; i < VtkVecSize; i++) {
				vtkVec[i] = (float) linalVec(i);
			}
			vtkArr->InsertNextTuple(vtkVec);
		};

		void insertQuantity(PhysicalQuantities::T quantity, vtkSmartPointer<vtkFloatArray> vtkArr,
		                           const int x, const int y, const int z) {
			auto Get = TGrid::PdeVector::QUANTITIES.at(quantity).Get;
			vtkArr->InsertNextValue((float)Get(grid->get(x, y, z)));
		};

		void insertOdeQuantity(PhysicalQuantities::T quantity, vtkSmartPointer<vtkFloatArray> vtkArr,
		                              const int x, const int y, const int z) {
			auto Get = TGrid::Model::InternalOde::QUANTITIES.at(quantity).Get;
			vtkArr->InsertNextValue((float)Get(grid->odeValues[grid->getIndex(x, y, z)]));
		};

		std::string makeFileNameForSnapshot(const int step);
	};
}

#endif // LIBGCM_VTKSTRUCTUREDSNAPSHOTTER_HPP
