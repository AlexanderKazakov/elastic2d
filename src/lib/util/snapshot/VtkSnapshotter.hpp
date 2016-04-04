#ifndef LIBGCM_VTKSNAPSHOTTER_HPP
#define LIBGCM_VTKSNAPSHOTTER_HPP

#include <lib/util/snapshot/Snapshotter.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>
#include <lib/mesh/grid/Cgal2DGrid.hpp>

#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkTriangle.h>

namespace gcm {
template<typename TMesh>
class VtkSnapshotter : public Snapshotter {
typedef typename TMesh::VtkGridType   VtkGridType;
typedef typename TMesh::VtkWriterType VtkWriterType;
const std::string FOLDER_NAME = std::string("vtk");

bool enableSnapshotting = false;
std::vector<PhysicalQuantities::T> quantitiesToSnap;

virtual void beforeStatementImpl(const Statement& statement) override {
	enableSnapshotting = statement.vtkSnapshotter.enableSnapshotting;
	quantitiesToSnap = statement.vtkSnapshotter.quantitiesToSnap;
}

virtual void snapshotImpl(const AbstractGrid* _mesh, const int step) override {
	if (!enableSnapshotting) { return; }
	mesh = dynamic_cast<const TMesh*>(_mesh);
	assert_true(mesh);
	// TODO - in this approach, we create a structure of size of the whole mesh,
	// and then write it, allocate size of mesh at every time step is not a good idea
	vtkGrid = VtkGridType::New();
	writeGeometry(*mesh, vtkGrid);

	for (auto& quantity : TMesh::PdeVector::VECTORS) {
		writeQuantity(quantity.first, &VtkSnapshotter::insertVector, 3);
	}
	for (auto& quantity : quantitiesToSnap) {
		writeQuantity(quantity, &VtkSnapshotter::insertQuantity, 1);
	}
	for (auto& quantity : TMesh::Model::InternalOde::QUANTITIES) {
		writeQuantity(quantity.first, &VtkSnapshotter::insertOdeQuantity, 1);
	}

	vtkWriter = VtkWriterType::New();
#ifdef CONFIG_VTK_5
	vtkWriter->SetInput(vtkGrid);
#else
	vtkWriter->SetInputData(vtkGrid);
#endif
	vtkWriter->SetFileName(makeFileNameForSnapshot(step,
	                                               vtkWriter->GetDefaultFileExtension(),
	                                               FOLDER_NAME).c_str());
	vtkWriter->Write();
}

USE_AND_INIT_LOGGER("gcm.VtkSnapshotter")
const TMesh * mesh;
vtkSmartPointer<VtkGridType> vtkGrid;
vtkSmartPointer<VtkWriterType> vtkWriter;

/**
 * The code below is for writing pde, ode, etc data for any rheology model
 */

typedef void (VtkSnapshotter::* InsertFunc)(PhysicalQuantities::T quantity,
                                            vtkSmartPointer<vtkFloatArray> vtkArr,
                                            const typename TMesh::VtkIterator& it);

void writeQuantity(PhysicalQuantities::T quantity, InsertFunc insertFunc,
                   const int numOfComponents) {
	auto vtkArr = vtkSmartPointer<vtkFloatArray>::New();
	vtkArr->SetNumberOfComponents(numOfComponents);
	vtkArr->Allocate((vtkIdType)mesh->sizeOfRealNodes(), 0);
	vtkArr->SetName(PhysicalQuantities::NAME.at(quantity).c_str());
	for (auto it = mesh->vtkBegin(); it != mesh->vtkEnd(); ++it) {
		(this->*insertFunc)(quantity, vtkArr, it);
	}
	vtkGrid->GetPointData()->AddArray(vtkArr);
}

void insertVector(PhysicalQuantities::T quantity, vtkSmartPointer<vtkFloatArray> vtkArr,
                  const typename TMesh::VtkIterator& it) {
	auto Get = TMesh::PdeVector::VECTORS.at(quantity).Get;
	auto linalVec = Get(mesh->pde(it));
	float vtkVec[3];
	for (int i = 0; i < 3; i++) {
		vtkVec[i] = (float) linalVec(i);
	}
	vtkArr->InsertNextTuple(vtkVec);
}

void insertQuantity(PhysicalQuantities::T quantity, vtkSmartPointer<vtkFloatArray> vtkArr,
                    const typename TMesh::VtkIterator& it) {
	auto Get = TMesh::PdeVector::QUANTITIES.at(quantity).Get;
	vtkArr->InsertNextValue((float)Get(mesh->pde(it)));
}

void insertOdeQuantity(PhysicalQuantities::T quantity, vtkSmartPointer<vtkFloatArray> vtkArr,
                       const typename TMesh::VtkIterator& it) {
	auto Get = TMesh::Model::InternalOde::QUANTITIES.at(quantity).Get;
	vtkArr->InsertNextValue((float)Get(mesh->ode(it)));
}

/**
 * Below are implementations of geometry writing for meshes of different types
 */
void writeGeometry(const CubicGrid& gcmGrid, vtkSmartPointer<vtkStructuredGrid> _vtkGrid) {
	auto sizes = gcmGrid.sizes;
	_vtkGrid->SetDimensions(sizes(0), sizes(1), sizes(2));

	auto points = vtkSmartPointer<vtkPoints>::New();
	points->Allocate((vtkIdType)gcmGrid.sizeOfRealNodes(), 0);
	for (auto it = gcmGrid.vtkBegin(); it != gcmGrid.vtkEnd(); ++it) {
		auto coords = gcmGrid.coords(it);
		real point[3] = {coords(0), coords(1), coords(2)};
		points->InsertNextPoint(point);
	}
	_vtkGrid->SetPoints(points);
}

void writeGeometry(const Cgal2DGrid& gcmGrid, vtkSmartPointer<vtkUnstructuredGrid> _vtkGrid) {
	auto points = vtkSmartPointer<vtkPoints>::New();
	points->Allocate((vtkIdType)gcmGrid.sizeOfRealNodes(), 0);
	for (auto it = gcmGrid.vtkBegin(); it != gcmGrid.vtkEnd(); ++it) {
		auto coords = gcmGrid.coords(it);
		real point[3] = {coords(0), coords(1), coords(2)};
		points->InsertNextPoint(point);
	}
	_vtkGrid->SetPoints(points);

	auto triangle = vtkSmartPointer<vtkTriangle>::New();
	for (auto it = gcmGrid.cellBegin(); it != gcmGrid.cellEnd(); ++it) {
		auto verticesIndices = gcmGrid.getVerticesOfCell(it);
		int i = 0;
		for (const auto& vI : verticesIndices) {
			triangle->GetPointIds()->SetId(i, (vtkIdType)vI);
			i++;
		}
		_vtkGrid->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds());
	}
}

};


}

#endif // LIBGCM_VTKSNAPSHOTTER_HPP
