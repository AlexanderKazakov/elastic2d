#ifndef LIBGCM_CGAL2DGRID_HPP
#define LIBGCM_CGAL2DGRID_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/lloyd_optimize_mesh_3.h>

#include <lib/mesh/UnstructuredGrid.hpp>
#include <lib/linal/linal.hpp>

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
		typedef CGAL::Delaunay_mesher_2<CDT, Criteria>              Mesher;
		typedef CDT::Vertex_handle                                  VertexHandle;
		typedef CDT::Face_handle                                    FaceHandle;
		typedef CDT::Geom_traits::Vector_2                          CgalVector2;
		typedef CDT::Point                                          Point;
		typedef CDT::Finite_faces_iterator                          FiniteFacesIterator;
		typedef CDT::Line_face_circulator                           LineFaceCirculator;
		
		typedef elements::Triangle<Iterator>                        Triangle;

		/**
		 * Iteration over all real mesh faces (there are also not meshing  ("!is_in_domain()") faces in the triangulation,
		 * this caused by the fact that CGAL triangulation is a convex hull of its vertices, we don't deal with them at all)
		 */
		struct RealFaceTester {
			bool operator()(const FiniteFacesIterator& it) const { return !it->is_in_domain(); };
		} realFaceTester;
		typedef CGAL::Filter_iterator<FiniteFacesIterator, RealFaceTester> CellIterator;
		CellIterator cellBegin() const { return CellIterator(triangulation.finite_faces_end(), realFaceTester, triangulation.finite_faces_begin()); };
		CellIterator cellEnd() const { return CellIterator(triangulation.finite_faces_end(), realFaceTester, triangulation.finite_faces_end()); };

	protected:
		/**
		 * @param it begin() <= iterator < end()
		 * @return index in std::vector
		 */
		size_t getIndex(const Iterator& it) const {
			return it.iter;
		};

	public:
		/** Read-only access to real coordinates with auxiliary 0 at z */
		const linal::Vector3 coords(const Iterator& it) const {
			auto point = vertexHandles[it.iter]->point();
			return {point.x(), point.y(), 0};
		};
		/** Read-only access to real coordinates */
		const linal::Vector2 coords2d(const Iterator& it) const {
			auto point = vertexHandles[it.iter]->point();
			return {point.x(), point.y()};
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

		virtual void initializeImpl(const Task &task) override;
		void triangulate();

		Triangle findOwnerTriangle(const Iterator& it, const linal::Vector2& shift) const {
			Triangle ans;
			auto ownerFace = findOwnerFace(it, CgalVector2(shift(0), shift(1)));
			if ( ownerFace->is_in_domain() && ownerFace != triangulation.infinite_face() ) {
				ans.inner = true;
				for (int i = 0; i < Triangle::N; i++) {
					ans.p[i] = Iterator((size_t)(verticesIndices.at(ownerFace->vertex(i))));
				}
			}
			return ans;
		};
		FaceHandle findOwnerFace(const Iterator& it, const CgalVector2 shift) const {
			auto beginVertex = vertexHandles[getIndex(it)];
			auto q = beginVertex->point() + shift; // point to find owner face for
			return triangulation.locate(q, beginVertex->incident_faces());
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
