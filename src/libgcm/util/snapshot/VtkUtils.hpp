#ifndef LIBGCM_VTKUTILS_HPP
#define LIBGCM_VTKUTILS_HPP

// disable warnings from vtk headers
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wlogical-op"


#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkTriangle.h>
#include <vtkTetra.h>
#include <vtkLine.h>

#pragma GCC diagnostic pop

#include <libgcm/util/Elements.hpp>


namespace gcm {

template<int> class CubicGrid;
template<int, template<int, typename, typename> class> class SimplexGrid;

namespace vtk_utils {

typedef elements::Element<Real2, 3> Triangle;
typedef elements::Element<Real3, 4> Tetra;

/// Types binding @{
template<template<int, typename, typename> class TriangulationT>
vtkSmartPointer<vtkTriangle> getVtkCell(const SimplexGrid<2, TriangulationT>&) {
	return vtkSmartPointer<vtkTriangle>::New();
}
template<template<int, typename, typename> class TriangulationT>
vtkSmartPointer<vtkTetra> getVtkCell(const SimplexGrid<3, TriangulationT>&) {
	return vtkSmartPointer<vtkTetra>::New();
}

template<int D>
vtkSmartPointer<vtkStructuredGrid> getVtkGrid(const CubicGrid<D>&) {
	return vtkSmartPointer<vtkStructuredGrid>::New();
}
template<int D, template<int, typename, typename> class TriangulationT>
vtkSmartPointer<vtkUnstructuredGrid> getVtkGrid(const SimplexGrid<D, TriangulationT>&) {
	return vtkSmartPointer<vtkUnstructuredGrid>::New();
}

inline vtkSmartPointer<vtkXMLStructuredGridWriter> getVtkWriter(
		vtkSmartPointer<vtkStructuredGrid>) {
	return vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
}
inline vtkSmartPointer<vtkXMLUnstructuredGridWriter> getVtkWriter(
		vtkSmartPointer<vtkUnstructuredGrid>) {
	return vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
}

template<typename VtkGridPtr>
std::string getVtkFileExtension(VtkGridPtr vtkGrid) {
	return getVtkWriter(vtkGrid)->GetDefaultFileExtension();
}
/// @}

/** Write cubic grid geometry to vtk */
template<int D>
void writeGeometry(const CubicGrid<D>& gcmGrid,
		vtkSmartPointer<vtkStructuredGrid> vtkGrid);

/** Write simplex grid geometry to vtk */
template<int D, template<int, typename, typename> class TriangulationT>
void writeGeometry(const SimplexGrid<D, TriangulationT>& gcmGrid,
		vtkSmartPointer<vtkUnstructuredGrid> vtkGrid);

/** Draw the line between the two points to vtk file */
void drawSegmentToVtk(
		const Real3 a, const Real3 b, const std::string& fileName = "line");

/** Draw the line between the two points to vtk file */
void drawSegmentToVtk(
		const Real2 a, const Real2 b, const std::string& fileName = "line");

/** Draw triangles to vtk file */
void drawCellsToVtk(
		const std::vector<Triangle>& cells, const std::string& fileName = "cells");

/** Draw tetrahedrons to vtk file */
void drawCellsToVtk(
		const std::vector<Tetra>& cells, const std::string& fileName = "cells");

/** Open the grid at runtime in vtk renderer (but Paraview is much more useful) */
void showGridInVtkRender(
		vtkSmartPointer<vtkUnstructuredGrid> grid);


template<typename VtkGridPtr>
void writeToFile(VtkGridPtr grid, const std::string& name) {
/// write vtk grid to file
	auto writer = getVtkWriter(grid);
#ifdef CONFIG_VTK_5
	writer->SetInput(grid);
#else
	writer->SetInputData(grid);
#endif
	writer->SetFileName(name.c_str());
	writer->Write();
}


template<typename GcmGrid>
void dumpGridToVtk(
/// write geometry of the given gcm grid to vtk file
		const GcmGrid& grid, const std::string& name = "grid") {
	auto vtkGrid = getVtkGrid(grid);
	writeGeometry(grid, vtkGrid);
	writeToFile(vtkGrid, name + "." + getVtkFileExtension(vtkGrid));
}


template<typename GcmGrid, typename VtkGrid>
void writeVertices(const GcmGrid& gcmGrid, vtkSmartPointer<VtkGrid> vtkGrid) {
/// write coordinates of the vertices of the grid
	auto points = vtkSmartPointer<vtkPoints>::New();
	points->Allocate((vtkIdType)gcmGrid.sizeOfRealNodes(), 0);
	for (auto it = gcmGrid.vtkBegin(); it != gcmGrid.vtkEnd(); ++it) {
		Real3 coords = gcmGrid.coords(it);
		real point[3] = {coords(0), coords(1), coords(2)};
		points->InsertNextPoint(point);
	}
	vtkGrid->SetPoints(points);
}


template<typename GcmGrid>
void writeCells(
		const GcmGrid& gcmGrid, vtkSmartPointer<vtkUnstructuredGrid> vtkGrid) {
/// write cells of the unstructured grid
	auto vtkCell = getVtkCell(gcmGrid);
	for (auto it = gcmGrid.cellBegin(); it != gcmGrid.cellEnd(); ++it) {
		const auto gcmCell = gcmGrid.createCell(*it);
		for (int i = 0; i < gcmCell.N; i++) {
			vtkCell->GetPointIds()->SetId(i, (vtkIdType)(gcmCell(i).iter));
		}
		vtkGrid->InsertNextCell(vtkCell->GetCellType(), vtkCell->GetPointIds());
	}
}


template<typename GcmMesh, typename VtkGridPtr, typename ValueInserter>
void addFieldToVertices(const int numberOfComponents, const std::string name, 
		const ValueInserter insertValue, const GcmMesh& gcmMesh, VtkGridPtr vtkGrid) {
/// write some field to the vertices of vtk grid from the vertices of gcm grid
	vtkSmartPointer<vtkFloatArray> vtkArr = vtkSmartPointer<vtkFloatArray>::New();
	vtkArr->SetNumberOfComponents(numberOfComponents);
	vtkArr->Allocate((vtkIdType)gcmMesh.sizeOfRealNodes());
	vtkArr->SetName(name.c_str());
	for (auto it = gcmMesh.vtkBegin(); it != gcmMesh.vtkEnd(); ++it) {
		insertValue(vtkArr, it);
	}
	vtkGrid->GetPointData()->AddArray(vtkArr);
}

} // namespace vtk_utils
} // namespace gcm

#endif // LIBGCM_VTKUTILS_HPP
