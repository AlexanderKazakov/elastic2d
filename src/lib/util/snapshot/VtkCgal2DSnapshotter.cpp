#include <lib/util/snapshot/VtkCgal2DSnapshotter.hpp>
#include <lib/grid/Cgal2DGrid.hpp>
#include <lib/rheology/models/Model.hpp>


using namespace gcm;

template<class TGrid>
void VtkCgal2DSnapshotter<TGrid>::snapshotImpl(const Grid* _grid, const int step) {
	LOG_DEBUG("Start snapshot writing to " << makeFileNameForSnapshot(step));
	grid = static_cast<const TGrid*>(_grid);

	vtkGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	grid->mesh.writeVtk(vtkGrid);

	/*vtkGrid->getPointData*/


	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
#ifdef CONFIG_VTK_5
	writer->SetInput(vtkGrid);
#else
	writer->SetInputData(vtkGrid);
#endif
	writer->SetFileName(makeFileNameForSnapshot(step).c_str());
	writer->Write();
}

template<class TGrid>
std::string VtkCgal2DSnapshotter<TGrid>::makeFileNameForSnapshot(const int step) {
	char buffer[50];
	sprintf(buffer, "%s%02d%s%05d%s", "snaps/core", MPI::COMM_WORLD.Get_rank(), "_snapshot", step, ".vtu");
	return std::string(buffer);
}


template class VtkCgal2DSnapshotter<Cgal2DGrid<Elastic2DModel>>;
template class VtkCgal2DSnapshotter<Cgal2DGrid<ContinualDamageElastic2DModel>>;
template class VtkCgal2DSnapshotter<Cgal2DGrid<IdealPlastic2DModel>>;
