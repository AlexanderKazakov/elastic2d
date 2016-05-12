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
	typedef CDT::Point                                          CgalPoint2;
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
	/** Iteration over all inner nodes */
	///@{
	typedef typename std::set<size_t>::const_iterator InnerIterator;
	InnerIterator innerBegin() const { return innerIndices.begin(); }
	InnerIterator innerEnd() const { return innerIndices.end(); }
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
		auto point = vertexHandles[getIndex(it)]->point();
		return {point.x(), point.y(), 0};
	}

	/** Read-only access to real coordinates */
	Real2 coords2d(const Iterator& it) const {
		return real2(vertexHandles[getIndex(it)]->point());
	}
	
	/** Border normal in specified point */
	Real2 normal(const Iterator& it) const {
		return normal(findBorderNeighbors(it));
	}
	
	/** All vertices connected with specified vertex */
	std::set<Iterator> findNeighborVertices(const Iterator& it) const;

	/**
	 * @param it begin() <= iterator < end()
	 * @return index in std::vector
	 */
	size_t getIndex(const Iterator& it) const {
		return it.iter;
	}

	/** @return indices of all vertices in vertexHandles which specified cell owns */
	std::array<size_t, 3> getVerticesOfCell(const CellIterator& it) const {
		std::array<size_t, 3> vertices;
		for (size_t i = 0; i < 3; i++) {
			vertices[i] = it->vertex((int)i)->info();
		}
		return vertices;
	}
	
	Iterator getIterator(const VertexHandle& v) const {
		return v->info();
	}
	
	bool isBorder(const Iterator& it) const {
		return borderIndices.find(getIndex(it)) != borderIndices.end();
	}

	size_t sizeOfRealNodes() const {
		return vertexHandles.size();
	}

	size_t sizeOfAllNodes() const {
		return sizeOfRealNodes();
	}

	real getMinimalSpatialStep() const {
		assert_gt(effectiveSpatialStep, 0);
		return effectiveSpatialStep;
	}

	/**
	 * Find triangle that contains point on specified distance (shift) 
	 * from specified point (it) by line walk from it to (it+shift).
	 * If the line goes out of the body, returned triangle.inner == false, 
	 * triangle points are not set.
	 * @note for convex bodies result is the same with locateOwnerTriangle
	 */
	Triangle findOwnerTriangle(const Iterator& it, const Real2& shift) const;
	
	/**
	 * Find triangle that contains point on specified distance (shift) 
	 * from specified point (it) by triangulation.locate function. 
	 * If found face is not "in_domain", returned triangle.inner == false, 
	 * triangle points are not set.
	 * @note for convex bodies result is the same with findOwnerTriangle
	 */
	Triangle locateOwnerTriangle(const Iterator& it, const Real2& shift) const;
	
	/**
	 * Starting from specified point along the line in specified direction,
	 * find the nearest border edge crossed by the line.
	 * Border edge represented by the pair of its vertices.
	 * @note cases when go from border node outside the body aren't handled
	 * @return crossed border edge as the pair of its border points
	 */
	std::pair<Iterator, Iterator> findCrossingBorder(
			const Iterator& start, const Real2& direction) const;
	
	/**
	 * @return pair of border vertices {v1, v2} incident to given border vertex;
	 * the order of border vertices is so that outer body normal = perpendicularClockwise(v1-v2)
	 * @see normal
	 */
	std::pair<Iterator, Iterator> findBorderNeighbors(const Iterator& it) const;
	
	/**
	 * Given with two neighbor border points, go along the border from first to 
	 * second and so on while the border is straight line collinear to (second - first)
	 * @return the point where the border bends from the line
	 */
	Iterator findBorderFlexion(Iterator first, Iterator second) const;


protected:
	/// Data
	///@{
	CDT triangulation;                       ///< CGAL triangulation data structure
	std::vector<VertexHandle> vertexHandles; ///< CGAL-"pointers" to each grid vertex
	std::set<size_t> borderIndices;          ///< indices of border vertices in vertexHandles
	std::set<size_t> innerIndices;           ///< indices of inner vertices in vertexHandles
	real effectiveSpatialStep = 0;           ///< used in triangulation criteria and Courant condition
	bool movable = false;                    ///< deformable(true) or immutable(false) grid
	///@}

	/** Move specified point on specified distance */
	void move(const Iterator& it, const Real2& d) {
		assert_true(movable);
		auto& point = vertexHandles[getIndex(it)]->point();
		point = point + cgalVector2(d);
	}

	
private:
	/** Functions for building the triangulation */
	///@{
	void triangulate(const Task::Cgal2DGrid& task);
	void insertPolygon(const Polygon& polygon);
	static Polygon makePolygon(const std::vector<Real2>& points);
	static CgalPoint2 findInnerPoint(const Polygon& polygon);
	void markInnersAndBorders();
	///@}

	/// Auxilliary functions for handling numerical methods queries 
	/// @{
	FaceHandle findCrossedFace(const VertexHandle start, const Real2& direction) const;
	
	std::pair<Iterator, Iterator> findCrossingBorder(
			const VertexHandle& v, const FaceHandle& f, const Real2& direction) const;
	
	VertexHandle commonVertex(const FaceHandle& a, const FaceHandle& b) const;

	Real2 normal(const std::pair<Iterator, Iterator>& borderNeighbors) const {
		VertexHandle first  = vertexHandles[getIndex(borderNeighbors.first)];
		VertexHandle second = vertexHandles[getIndex(borderNeighbors.second)];
		Real2 borderVector = real2(first->point()) - real2(second->point());
		return linal::normalize(
		       linal::perpendicularClockwise(borderVector));
	}
	/// @}


	void printFace(const FaceHandle& f, const std::string& name = "") const;
 
	static CgalPoint2 cgalPoint2(const Real2& p) {
		return CgalPoint2(p(0), p(1));
	}
	
	static Real2 real2(const CgalPoint2& p) {
		return {p.x(), p.y()};
	}
	
	static CgalVector2 cgalVector2(const Real2& p) {
		return CgalVector2(p(0), p(1));
	}

	USE_AND_INIT_LOGGER("gcm.Cgal2DGrid")
};


}

#endif // LIBGCM_CGAL2DGRID_HPP
