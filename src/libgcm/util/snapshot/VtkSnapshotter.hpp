#ifndef LIBGCM_VTKSNAPSHOTTER_HPP
#define LIBGCM_VTKSNAPSHOTTER_HPP

#include <libgcm/util/snapshot/Snapshotter.hpp>
#include <libgcm/util/snapshot/VtkUtils.hpp>

namespace gcm {

template<typename TMesh>
class VtkSnapshotter : public Snapshotter {
public:
	typedef typename TMesh::Grid           Grid;
	typedef typename TMesh::VtkIterator    VtkIter;
	const std::string FOLDER_NAME = std::string("vtk");
	
	VtkSnapshotter(const Task& task) :
			Snapshotter(task),
			quantitiesToSnap(task.vtkSnapshotter.quantitiesToSnap) { }
	
	virtual void snapshotImpl(const AbstractGrid* _mesh, const int step) override {
		const TMesh* mesh = dynamic_cast<const TMesh*>(_mesh);
		assert_true(mesh);
		// note: we allocate a structure of size of the whole mesh at every time step!
		auto vtkGrid = vtk_utils::getVtkGrid(*mesh);
		vtk_utils::writeGeometry(*mesh, vtkGrid);
	
		for (auto quantity : TMesh::Model::PdeVariables::VECTORS) {
			auto Get = TMesh::Model::PdeVariables::VECTORS.at(quantity.first).Get;
			vtk_utils::addFieldToVertices(
					3, PhysicalQuantities::NAME.at(quantity.first),
					[&](vtkSmartPointer<vtkFloatArray> vtkArr, VtkIter it) {
						auto linalVec = Get(mesh->pde(it));
						float vtkVec[3];
						for (int i = 0; i < 3; i++) {
							vtkVec[i] = (float)linalVec(i);
						}
						vtkArr->InsertNextTuple(vtkVec);
					},
					*mesh, vtkGrid);
		}
		for (PhysicalQuantities::T quantity : quantitiesToSnap) {
			auto Get = TMesh::Model::PdeVariables::QUANTITIES.at(quantity).Get;
			vtk_utils::addFieldToVertices(
					1, PhysicalQuantities::NAME.at(quantity),
					[&](vtkSmartPointer<vtkFloatArray> vtkArr, VtkIter it) {
						vtkArr->InsertNextValue((float)Get(mesh->pde(it)));
					},
					*mesh, vtkGrid);
		}
		vtk_utils::addFieldToVertices(
				1, "material_index",
				[&](vtkSmartPointer<vtkFloatArray> vtkArr, VtkIter it) {
					vtkArr->InsertNextValue(
							(float)(mesh->material(it)->materialNumber));
				},
				*mesh, vtkGrid);
		
		vtk_utils::writeToFile(vtkGrid, makeFileNameForSnapshot(
				std::to_string(mesh->id), step,
				vtk_utils::getVtkFileExtension(vtkGrid), FOLDER_NAME));
	}
	
	
private:
	const std::vector<PhysicalQuantities::T> quantitiesToSnap;
	USE_AND_INIT_LOGGER("gcm.VtkSnapshotter")
	
};


}

#endif // LIBGCM_VTKSNAPSHOTTER_HPP
