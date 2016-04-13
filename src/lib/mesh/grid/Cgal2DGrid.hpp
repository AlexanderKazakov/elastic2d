#ifndef LIBGCM_CGAL2DGRID_HPP
#define LIBGCM_CGAL2DGRID_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Polygon_2.h>
//#include <CGAL/Polygon_2_algorithms.h>

#include <lib/mesh/grid/UnstructuredGrid.hpp>
#include <lib/linal/linal.hpp>

namespace gcm {
/**
 * 2D movable unstructured triangle grid by CGAL library
 */
class Cgal2DGrid : public UnstructuredGrid {
public:
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Polygon_2<K>                                  Polygon;
	typedef CGAL::Triangulation_vertex_base_with_info_2
	                                                <size_t, K> Vb;
	typedef CGAL::Delaunay_mesh_face_base_2<K>                  Fb;
	typedef CGAL::Triangulation_data_structure_2<Vb, Fb>        Tds;
	typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>  CDT;
	typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>            Criteria;
	typedef CGAL::Delaunay_mesher_2<CDT, Criteria>              Mesher;
	typedef CDT::Vertex_handle                                  VertexHandle;
	typedef CDT::Face_handle                                    FaceHandle;
	typedef CDT::Geom_traits::Vector_2                          CgalVector2;
	typedef CDT::Triangle                                       CgalTriangle2;
	typedef CDT::Point                                          Point;
	typedef CDT::Finite_faces_iterator                          FiniteFacesIterator;
	typedef CDT::Finite_vertices_iterator                       FiniteVerticesIterator;
	typedef CDT::Line_face_circulator                           LineFaceCirculator;

	typedef elements::Triangle<Iterator>                        Triangle;
	
	/// @name Iterators 
	///@{

	typedef Iterator ForwardIterator;
	ForwardIterator begin() const { return 0; }
	ForwardIterator end() const { return sizeOfRealNodes(); }
	typedef Iterator VtkIterator;
	VtkIterator vtkBegin() const { return begin(); }
	VtkIterator vtkEnd() const { return end(); }
	
	/** Iteration over all border nodes */
	///@{
	typedef typename std::set<size_t>::const_iterator BorderIterator;
	
	BorderIterator borderBegin() const { return borderIndices.begin(); }
	BorderIterator borderEnd() const { return borderIndices.end(); }
	///@}

	/**
	 * Iteration over all real mesh faces.
	 * "Real" means face which covers the body not outer space ("is_in_domain()").
	 * There are also not meshing ("!is_in_domain()") faces in the triangulation.
	 * They cover the space around bodies, because CGAL triangulation 
	 * is a convex hull of its vertices.
	 */
	///@{
	struct RealFaceTester {
		bool operator()(const FiniteFacesIterator& it) const {
			return !it->is_in_domain();
		}
	} realFaceTester;
	
	typedef CGAL::Filter_iterator<FiniteFacesIterator, RealFaceTester> CellIterator;
	
	CellIterator cellBegin() const { 
		return CellIterator(triangulation.finite_faces_end(),
		                    realFaceTester,
		                    triangulation.finite_faces_begin()); }
	CellIterator cellEnd() const { 
		return CellIterator(triangulation.finite_faces_end(),
		                    realFaceTester,
		                    triangulation.finite_faces_end()); }
	///@}
	///@}

	Cgal2DGrid(const Task& task);
	virtual ~Cgal2DGrid() { }

	/** Read-only access to real coordinates with auxiliary 0 at z */
	Real3 coords(const Iterator& it) const {
		auto point = vertexHandles[it.iter]->point();
		return {point.x(), point.y(), 0};
	}

	/** Read-only access to real coordinates */
	Real2 coords2d(const Iterator& it) const {
		return PointToReal2(vertexHandles[it.iter]->point());
	}
	
	/** Border normal in specified point */
	Real2 normal(const Iterator& iter) const;

protected:
	/** Move specified point on specified distance */
	void move(const Iterator& it, const Real2& d) {
		assert_true(movable);
		auto& point = vertexHandles[it.iter]->point();
		point = point + CgalVector2(d(0), d(1));
	}

	/**
	 * @param it begin() <= iterator < end()
	 * @return index in std::vector
	 */
	size_t getIndex(const Iterator& it) const {
		return it.iter;
	}

public:
	/** @return indices of all vertices in vertexHandles which specified cell owns */
	std::array<size_t, 3> getVerticesOfCell(const CellIterator& it) const {
		std::array<size_t, 3> vertices;
		for (size_t i = 0; i < 3; i++) {
			vertices[i] = it->vertex((int)i)->info();
		}
		return vertices;
	}
	
	const CDT& getTriangulation() const {
		return triangulation;
	}
	
	Iterator getIterator(const VertexHandle& v) const {
		return v->info();
	}

	size_t sizeOfRealNodes() const {
		return vertexHandles.size();
	}

	size_t sizeOfAllNodes() const {
		return sizeOfRealNodes();
	}

protected:
	/// Data
	///@{
	CDT triangulation;                               ///< CGAL triangulation data structure
	std::vector<VertexHandle> vertexHandles;         ///< CGAL-"pointers" to each grid vertex
	std::set<size_t> borderIndices;                  ///< indices of border vertices in vertexHandles
	real effectiveSpatialStep = 0;                   ///< used in triangulation criteria and Courant condition
	bool movable = false;                            ///< deformable(true) or immutable(false) grid
	///@}

	real getMinimalSpatialStep() const {
		assert_gt(effectiveSpatialStep, 0);
		return effectiveSpatialStep;
	}

	Triangle findOwnerTriangle(const Iterator& it, const Real2& shift) const {
		Triangle ans;
		auto ownerFace = findOwnerFace(it, CgalVector2(shift(0), shift(1)));
		if (ownerFace->is_in_domain()) {
			ans.inner = true;
			for (int i = 0; i < Triangle::N; i++) {
				ans.p[i] = getIterator(ownerFace->vertex(i));
			}
		}
		return ans;
	}

	FaceHandle findOwnerFace(const Iterator& it, const CgalVector2 shift) const {
		auto beginVertex = vertexHandles[getIndex(it)];
		auto q = beginVertex->point() + shift; // point to find owner face for
		return triangulation.locate(q, beginVertex->incident_faces());
	}

private:
	/** Functions for building the triangulation */
	///@{
	void triangulate(const Task::Cgal2DGrid& task);
	void insertPolygon(const Polygon& polygon);
	static Polygon makePolygon(const std::vector<Real2>& points);
	static Point findInnerPoint(const Polygon& polygon);
	void markBorders();
	///@}
 
	static Point Real2ToPoint(const Real2& p) {
		return Point(p(0), p(1));
	}
	
	static Real2 PointToReal2(const Point& p) {
		return {p.x(), p.y()};
	}

	USE_AND_INIT_LOGGER("gcm.Cgal2DGrid")
};


}

#endif // LIBGCM_CGAL2DGRID_HPP
