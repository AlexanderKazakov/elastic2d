#ifndef LIBGCM_CGAL2DGRID_HPP
#define LIBGCM_CGAL2DGRID_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <lib/mesh/UnstructuredGrid.hpp>
#include <lib/linal/linal.hpp>
#include <lib/util/areas/AxisAlignedBoxArea.hpp>

namespace gcm {
	/**
	 * 2D movable unstructured triangle grid by CGAL library
	 */
	class Cgal2DGrid : public UnstructuredGrid {
	public:
		virtual ~Cgal2DGrid() { };
		typedef Iterator ForwardIterator;
		ForwardIterator begin() const { return 0; };
		ForwardIterator end() const { return sizeOfRealNodes(); };
		typedef Iterator VtkIterator;
		VtkIterator vtkBegin() const { return begin(); };
		VtkIterator vtkEnd() const { return end(); };

		typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
		typedef CGAL::Triangulation_vertex_base_2<K>                Vb;
		typedef CGAL::Delaunay_mesh_face_base_2<K>                  Fb;
		typedef CGAL::Triangulation_data_structure_2<Vb, Fb>        Tds;
		typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>  CDT;
		typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>            Criteria;
		typedef CDT::Vertex_handle                                  VertexHandle;
		typedef CDT::Face_handle                                    FaceHandle;
		typedef CDT::Point                                          Point;
		typedef CDT::Finite_faces_iterator                          CellIterator;

		CellIterator cellBegin() const { return triangulation.finite_faces_begin(); };
		CellIterator cellEnd() const { return triangulation.finite_faces_end(); };

	protected:
		/**
		 * @param it begin() <= iterator < end()
		 * @return index in std::vector
		 */
		size_t getIndex(const Iterator& it) const {
			return it.iter;
		};

	public:
		/** Read-only access to real coordinates */
		const linal::Vector3 coords(const Iterator& it) const {
			auto point = vertexHandles[it.iter]->point();
			return {point.x(), point.y(), 0};
		};

		/** @return indices of all vertices in vertexHandles which specified cell owns */
		void getVerticesOfCell(const CellIterator& it, int (&vertices)[3]) const {
			for( int i = 0; i < 3; i++) {
				vertices[i] = verticesIndices.at(it->vertex(i));
			}
		};

		size_t sizeOfRealNodes() const {
			return vertexHandles.size();
		};
		size_t sizeOfAllNodes() const {
			return sizeOfRealNodes();
		};

	protected:
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

		virtual void initializeImpl(const Task &task) override {
			LOG_INFO("Start initialization");
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

			recalculateMinimalSpatialStep();
			initializeImplImpl(task);
		};

		virtual void beforeStageImpl() override { };
		virtual void afterStageImpl() override { };
		virtual void beforeStepImpl() override { };
		virtual void afterStepImpl() override { };
		virtual void recalculateMinimalSpatialStep() override { minimalSpatialStep = 1; /* TODO  */ };

		virtual void applyInitialConditions(const Task& task) = 0;
		virtual void recalculateMaximalLambda() = 0;

		virtual void initializeImplImpl(const Task& task) = 0;

		USE_AND_INIT_LOGGER("gcm.Cgal2DGrid");
	};
}

#endif // LIBGCM_CGAL2DGRID_HPP
