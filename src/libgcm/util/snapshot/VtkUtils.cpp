#include <libgcm/util/snapshot/VtkUtils.hpp>
#include <libgcm/grid/simplex/cgal/CgalTriangulation.hpp>
#include <libgcm/grid/cubic/CubicGrid.hpp>
#include <libgcm/grid/simplex/SimplexGrid.hpp>

// disable warnings from vtk headers
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"

#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

#pragma GCC diagnostic pop


namespace gcm {
namespace vtk_utils {

template<int D>
void writeGeometry(const CubicGrid<D>& gcmGrid,
		vtkSmartPointer<vtkStructuredGrid> vtkGrid) {
	Int3 sizes = Int3::Ones();
	for (int i = 0; i < D; i++) {
		sizes(i) = gcmGrid.sizes(i);
	}
	vtkGrid->SetDimensions(sizes(0), sizes(1), sizes(2));
	writeVertices(gcmGrid, vtkGrid);
}


template<int D, template<int, typename, typename> class TriangulationT>
void writeGeometry(const SimplexGrid<D, TriangulationT>& gcmGrid,
		vtkSmartPointer<vtkUnstructuredGrid> vtkGrid) {
	writeVertices(gcmGrid, vtkGrid);
	writeCells(gcmGrid, vtkGrid);
	typedef typename SimplexGrid<D, CgalTriangulation>::VtkIterator VtkIter;
	vtk_utils::addFieldToVertices(
			1, "index_of_node",
			[&](vtkSmartPointer<vtkFloatArray> vtkArr, VtkIter it) {
				vtkArr->InsertNextValue((float)(gcmGrid.getIndex(it)));
			},
			gcmGrid, vtkGrid);
}


void drawCellsToVtk(
		const std::vector<Tetra>& cells, const std::string& fileName) {
	const std::vector<Real3> pointsVec =
			elements::sortedUniquePointsOfElements(cells);
	vtkSmartPointer<vtkPoints> pointsVtk = vtkSmartPointer<vtkPoints>::New();
	for (const Real3& point : pointsVec) {
		pointsVtk->InsertNextPoint(point(0), point(1), point(2));
	}
	
	vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
	for (const Tetra& cell : cells) {
		vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();
		assert_eq(cell.N, cell.n);
		for (int i = 0; i < cell.n; i++) {
			tetra->GetPointIds()->SetId(i,
					(vtkIdType)Utils::findIndexOfValueInSortedArray(
							pointsVec.begin(), pointsVec.end(), cell(i)));
		}
		cellArray->InsertNextCell(tetra);
	}
	
	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	grid->SetPoints(pointsVtk);
	grid->SetCells(VTK_TETRA, cellArray);
	writeToFile(grid, fileName + "." + getVtkFileExtension(grid));
}


void drawCellsToVtk(
		const std::vector<Triangle>& cells, const std::string& fileName) {
	const std::vector<Real2> pointsVec =
			elements::sortedUniquePointsOfElements(cells);
	vtkSmartPointer<vtkPoints> pointsVtk = vtkSmartPointer<vtkPoints>::New();
	for (const Real2& point : pointsVec) {
		pointsVtk->InsertNextPoint(point(0), point(1), 0);
	}
	
	vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
	for (const Triangle& cell : cells) {
		vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
		assert_eq(cell.N, cell.n);
		for (int i = 0; i < cell.n; i++) {
			triangle->GetPointIds()->SetId(i,
					(vtkIdType)Utils::findIndexOfValueInSortedArray(
							pointsVec.begin(), pointsVec.end(), cell(i)));
		}
		cellArray->InsertNextCell(triangle);
	}
	
	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	grid->SetPoints(pointsVtk);
	grid->SetCells(VTK_TRIANGLE, cellArray);
	writeToFile(grid, fileName + "." + getVtkFileExtension(grid));
}


void drawSegmentToVtk(const Real3 a, const Real3 b, const std::string& fileName) {
	vtkSmartPointer<vtkPoints> pointsVtk = vtkSmartPointer<vtkPoints>::New();
	pointsVtk->InsertNextPoint(a(0), a(1), a(2));
	pointsVtk->InsertNextPoint(b(0), b(1), b(2));
	
	vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
	line->GetPointIds()->SetId(0, 0);
	line->GetPointIds()->SetId(1, 1);
	cellArray->InsertNextCell(line);
	
	vtkSmartPointer<vtkUnstructuredGrid> grid =
			vtkSmartPointer<vtkUnstructuredGrid>::New();
	grid->SetPoints(pointsVtk);
	grid->SetCells(VTK_LINE, cellArray);
	
	writeToFile(grid, fileName + "." + getVtkFileExtension(grid));
}


void drawSegmentToVtk(const Real2 a, const Real2 b, const std::string& fileName) {
	Real3 a_ = {a(0), a(1), 0}; Real3 b_ = {b(0), b(1), 0};
	drawSegmentToVtk(a_, b_, fileName);
}


void showGridInVtkRender(vtkSmartPointer<vtkUnstructuredGrid> grid) {
	// Create a mapper and actor
	vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(grid->GetProducerPort());
#else
	mapper->SetInputData(grid);
#endif
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	
	// Create a renderer, render window, and interactor
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = 
			vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	
	// Add the actor to the scene
	renderer->AddActor(actor);
	renderer->SetBackground(0.5, 0.5, 0.5); // Background color grey
	
	// Render and interact
	renderWindow->Render();
	renderWindowInteractor->Start();
}


template void writeGeometry(const CubicGrid<1>& gcmGrid,
		vtkSmartPointer<vtkStructuredGrid> vtkGrid);
template void writeGeometry(const CubicGrid<2>& gcmGrid,
		vtkSmartPointer<vtkStructuredGrid> vtkGrid);
template void writeGeometry(const CubicGrid<3>& gcmGrid,
		vtkSmartPointer<vtkStructuredGrid> vtkGrid);
template void writeGeometry(const SimplexGrid<2, CgalTriangulation>& gcmGrid,
		vtkSmartPointer<vtkUnstructuredGrid> vtkGrid);
template void writeGeometry(const SimplexGrid<3, CgalTriangulation>& gcmGrid,
		vtkSmartPointer<vtkUnstructuredGrid> vtkGrid);

} // namespace vtk_utils
} // namespace gcm
