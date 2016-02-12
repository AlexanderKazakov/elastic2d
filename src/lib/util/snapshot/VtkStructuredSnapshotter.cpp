#include <lib/util/snapshot/VtkStructuredSnapshotter.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/rheology/models/Model.hpp>


using namespace gcm;

template<class TGrid>
void VtkStructuredSnapshotter<TGrid>::snapshotImpl(const Grid* _grid, const int step) {
	LOG_DEBUG("Start snapshot writing to " << makeFileNameForSnapshot(step));
	grid = static_cast<const TGrid*>(_grid);

	vtkStrGrid = vtkSmartPointer<vtkStructuredGrid>::New();
	vtkStrGrid->SetDimensions(grid->sizes(0), grid->sizes(1), grid->sizes(2));

	// TODO - in this approach, we create a structure of size of the whole mesh,
	// and then write it, allocate size of mesh at every time step is not a good idea
	auto points = vtkSmartPointer<vtkPoints>::New();
	points->Allocate(linal::directProduct(grid->sizes), 0);
	for (int z = 0; z < grid->sizes(2); z++) {
		for (int y = 0; y < grid->sizes(1); y++) {
			for (int x = 0; x < grid->sizes(0); x++) {
				auto coords = grid->getCoordinates(x, y, z);
				real point[3] = {coords(0), coords(1), coords(2)};
				points->InsertNextPoint(point);
			}
		}
	}
	vtkStrGrid->SetPoints(points);

	for (auto& quantity : TGrid::PdeVector::VECTORS) {
		writeQuantity(quantity.first, &VtkStructuredSnapshotter::insertVector, VtkVecSize);
	}
	for (auto& quantity : quantitiesToWrite) {
		writeQuantity(quantity, &VtkStructuredSnapshotter::insertQuantity, 1);
	}
	for (auto& quantity : TGrid::Model::InternalOde::QUANTITIES) {
		writeQuantity(quantity.first, &VtkStructuredSnapshotter::insertOdeQuantity, 1);
	}

	auto writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
#ifdef CONFIG_VTK_5
	writer->SetInput(vtkStrGrid);
#else
	writer->SetInputData(vtkStrGrid);
#endif
	writer->SetFileName(makeFileNameForSnapshot(step).c_str());
	writer->Write();
}

template<class TGrid>
std::string VtkStructuredSnapshotter<TGrid>::makeFileNameForSnapshot(const int step) {
	char buffer[50];
	sprintf(buffer, "%s%02d%s%05d%s", "snaps/core", MPI::COMM_WORLD.Get_rank(), "_snapshot", step, ".vts");
	return std::string(buffer);
}


template class VtkStructuredSnapshotter<StructuredGrid<Elastic1DModel>>;
template class VtkStructuredSnapshotter<StructuredGrid<Elastic2DModel>>;
template class VtkStructuredSnapshotter<StructuredGrid<Elastic3DModel>>;
template class VtkStructuredSnapshotter<StructuredGrid<OrthotropicElastic3DModel>>;
template class VtkStructuredSnapshotter<StructuredGrid<ContinualDamageElastic2DModel>>;
template class VtkStructuredSnapshotter<StructuredGrid<IdealPlastic2DModel>>;
template class VtkStructuredSnapshotter<StructuredGrid<SuperDuperModel>>;
