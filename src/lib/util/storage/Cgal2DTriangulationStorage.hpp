#ifndef LIBGCM_CGAL2DTRIANGULATIONSTORAGE_HPP
#define LIBGCM_CGAL2DTRIANGULATIONSTORAGE_HPP

#include <vector>

#include <lib/util/areas/AxisAlignedBoxArea.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkTriangle.h>

namespace gcm {
	struct Cgal2DTriangulationStorage {

		typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
		typedef CGAL::Triangulation_vertex_base_2<K>                Vb;
		typedef CGAL::Delaunay_mesh_face_base_2<K>                  Fb;
		typedef CGAL::Triangulation_data_structure_2<Vb, Fb>        Tds;
		typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>  CDT;
		typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>            Criteria;
		typedef CDT::Vertex_handle                                  VertexHandle;
		typedef CDT::Face_handle                                    FaceHandle;
		typedef CDT::Point                                          Point;


		void initialize(const Task& task) {
			auto h = linal::plainDivision(task.lengthes, task.sizes - linal::VectorInt<3>({1, 1, 1}));
			auto sizes = task.sizes;
			auto startR = task.startR;
			auto finishR = startR + linal::plainMultiply(sizes, h);
			AxisAlignedBoxArea box(startR, finishR);
			triangulate(box, h);
			vertexHandles.resize(triangulation.number_of_vertices());
			size_t vertexIndex = 0;
			for (auto it = triangulation.finite_vertices_begin(); it != triangulation.finite_vertices_end(); it++) {
				auto handle = it->handle();
				vertexHandles[vertexIndex] = handle;
				verticesIndices.insert({handle, vertexIndex});
				vertexIndex++;
			}
		};

		linal::Vector3 getCoordinates(const size_t i) const {
			auto point = vertexHandles[i]->point();
			return {point.x(), point.y(), 0};
		};

		real getMinimalSpatialStep() const {
			return 1; // TODO
		};

		size_t verticesSize() const {
			return vertexHandles.size();
		};

		void writeVtk(vtkSmartPointer<vtkUnstructuredGrid> vtkGrid) const
		{
			auto points = vtkSmartPointer<vtkPoints>::New();
			points->Allocate((vtkIdType) verticesSize(), 0);
			for (size_t i = 0; i < verticesSize(); i++) {
				auto coords = getCoordinates(i);
				real point[3] = {coords(0), coords(1), coords(2)};
				points->InsertNextPoint(point);
				std::cout << coords;
			}
			vtkGrid->SetPoints(points);

			auto triangle = vtkSmartPointer<vtkTriangle>::New();
			for (auto it = triangulation.finite_faces_begin(); it != triangulation.finite_faces_end(); it++) {
				for( int i = 0; i < 3; i++) {
					int vertexIndex = verticesIndices.at(it->vertex(i));
					triangle->GetPointIds()->SetId(i, vertexIndex);
				}
				vtkGrid->InsertNextCell(triangle->GetCellType(),triangle->GetPointIds());
			}
		};

	private:
		CDT                         triangulation;
		std::vector<VertexHandle>   vertexHandles;
		std::map<VertexHandle, int> verticesIndices;

		void triangulate(const AxisAlignedBoxArea& box, linal::Vector<3> h) {
			auto startR = box.getMin();
			auto finishR = box.getMax();
			VertexHandle va = triangulation.insert(Point(startR(0), startR(1)));
			VertexHandle vb = triangulation.insert(Point(startR(0), finishR(1)));
			VertexHandle vc = triangulation.insert(Point(finishR(0), finishR(1)));
			VertexHandle vd = triangulation.insert(Point(finishR(0), startR(1)));
			triangulation.insert_constraint(va, vb);
			triangulation.insert_constraint(vb, vc);
			triangulation.insert_constraint(vc, vd);
			triangulation.insert_constraint(vd, va);
			std::cout << "Number of vertices: " << triangulation.number_of_vertices() << std::endl;
			std::cout << triangulation.number_of_faces() << std::endl;
			std::cout << "Meshing the triangulation..." << std::endl;
			CGAL::refine_Delaunay_mesh_2(triangulation, Criteria(0.125, h(0)));
			std::cout << "Number of vertices: " << triangulation.number_of_vertices() << std::endl;
			std::cout << triangulation.number_of_faces() << std::endl;
		};
	};
}

#endif // LIBGCM_CGAL2DTRIANGULATIONSTORAGE_HPP
