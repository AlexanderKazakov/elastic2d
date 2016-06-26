#ifndef LIBGCM_CGAL2DGRID_HPP
#define LIBGCM_CGAL2DGRID_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>


#include <lib/mesh/grid/UnstructuredGrid.hpp>

namespace CGAL {
	template <class Traits_P, class Container_P> class Polygon_2;
	template <class CDT> class Delaunay_mesh_size_criteria_2;
	template <typename Tr, typename Crit> class Delaunay_mesher_2;
}

namespace gcm {
/**
 * 2D movable unstructured triangle grid by CGAL library
 */
class Cgal2DGrid : public UnstructuredGrid {
public:
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
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
	typedef CGAL::Polygon_2<K, std::vector<CgalPoint2>>         Polygon;
	typedef CDT::Finite_faces_iterator                          FiniteFacesIterator;
	typedef CDT::Finite_vertices_iterator                       FiniteVerticesIterator;
	typedef CDT::Line_face_circulator                           LineFaceCirculator;

	typedef elements::Triangle<Iterator>                        Cell;
	
	/// Space dimensionality
	static const int DIMENSIONALITY = 2;
	/// An *estimation* of maximal possible number of vertices connected 
	/// with some inner vertex (it can be more in a very rare cases)
	static const int MAX_NUMBER_OF_NEIGHBOR_VERTICES = 8;
	
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
	typedef typename std::vector<size_t>::const_iterator BorderIterator;
	BorderIterator borderBegin() const { return borderIndices.begin(); }
	BorderIterator borderEnd() const { return borderIndices.end(); }
	///@}
	/** Iteration over all inner nodes */
	///@{
	typedef typename std::vector<size_t>::const_iterator InnerIterator;
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
		auto point = vertexHandle(it)->point();
		return {point.x(), point.y(), 0};
	}

	/** Read-only access to real coordinates */
	Real2 coordsD(const Iterator& it) const {
		return real2(vertexHandle(it)->point());
	}
	
	/** Border normal in specified point */
	Real2 normal(const Iterator& it) const {
		return normal(it, findBorderNeighbors(it));
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
		VertexHandle vh = vertexHandle(it);
		auto beginFace = triangulation.incident_faces(vh);
		auto faceCirculator = beginFace;
		do {
			if ( !faceCirculator->is_in_domain() ) {
				return true;
			}
			++faceCirculator;
		} while (faceCirculator != beginFace);
		
		return false;
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
	 * In some bad degenerate cases cell.n == 0 and points are not set.
	 * If the line goes out of the body immediately (from border in outer 
	 * direction), cell.n == 0 and points are not set.
	 * If the line goes through the body, but go out before reaching
	 * the target point, returned cell.n == 2, and only 2 points are set
	 * - they are two border vertices - the edge crossed by the line.
	 * If the line goes through the body, and reach the target point,
	 * returned cell.n == 3, and 3 points are set.
	 * @note for convex bodies and (n == 3) result is the same with locateOwnerCell
	 */
	Cell findOwnerCell(const Iterator& it, const Real2& shift) const;
	
	/**
	 * Locate triangle that contains point on specified distance (shift) 
	 * from specified point (it) by triangulation.locate function. 
	 * If found face is not "in_domain" (out of body), returned cell.n == 0, 
	 * and points are not set. Else n == 3 and all 3 points are set.
	 * @note for convex bodies and (n == 3) result is the same with findOwnerCell
	 */
	Cell locateOwnerCell(const Iterator& it, const Real2& shift) const;
	
	/**
	 * @return pair of border vertices {v1, v2} incident to given border vertex;
	 * the order of border vertices is so that outer body normal is met
	 * when going around given border vertex clockwise from v1 to v2
	 * @see normal
	 */
	std::pair<Iterator, Iterator> findBorderNeighbors(const Iterator& it) const;
		
	/** Find vertex with specified coordinates. Throw Exception if there isn't such. */
	Iterator findVertexByCoordinates(const Real2& coordinates) const;


protected:
	/// Data
	///@{
	CDT triangulation;                       ///< CGAL triangulation data structure
	std::vector<VertexHandle> vertexHandles; ///< CGAL-"pointers" to each grid vertex
	std::vector<size_t> borderIndices;       ///< indices of border vertices in vertexHandles
	std::vector<size_t> innerIndices;        ///< indices of inner vertices in vertexHandles
	real effectiveSpatialStep = 0;           ///< used in triangulation criteria and Courant condition
	bool movable = false;                    ///< deformable(true) or immutable(false) grid
	///@}

	/** Move specified point on specified distance */
	void move(const Iterator& it, const Real2& d) {
		assert_true(movable);
		auto& point = vertexHandle(it)->point();
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
	VertexHandle vertexHandle(const Iterator& it) const {
		return vertexHandles [getIndex(it)];
	}
	
	std::vector<VertexHandle> commonVertices(const FaceHandle& a, const FaceHandle& b) const;	

	Real2 normal(const Iterator& it,
			const std::pair<Iterator, Iterator>& borderNeighbors) const {
		VertexHandle v = vertexHandle(it);
		VertexHandle first  = vertexHandle(borderNeighbors.first);
		VertexHandle second = vertexHandle(borderNeighbors.second);
		
		Real2 borderVector1 = real2(first->point()) - real2(v->point());
		Real2 borderVector2 = real2(v->point()) - real2(second->point());
		
		Real2 norm1 = linal::normalize(
		              linal::perpendicularClockwise(borderVector1));
		Real2 norm2 = linal::normalize(
		              linal::perpendicularClockwise(borderVector2));
		
		return linal::normalize(norm1 + norm2);
	}
	
	bool contains(const FaceHandle face, const Real2& q) const {
	/// is triangle with small layer around contains the point 
		Real2 a = real2(face->vertex(0)->point());
		Real2 b = real2(face->vertex(1)->point());
		Real2 c = real2(face->vertex(2)->point());
		Real3 lambda = linal::barycentricCoordinates(a, b, c, q);
		return (lambda(0) > -EQUALITY_TOLERANCE) &&
		       (lambda(1) > -EQUALITY_TOLERANCE) &&
		       (lambda(2) > -EQUALITY_TOLERANCE);
	}
	
	Cell createCell(const Iterator& it, 
			const FaceHandle current, const FaceHandle previous) const {
	/// create Cell used as answer to numerical method queries about point location
		Cell ans;
		ans.n = 0;
		if (current == NULL) {
		// not found
			assert_true(previous == NULL);
			return ans;
		}
		
		if ( current->is_in_domain() ) {
			// return found cell
			assert_true(previous->is_in_domain());
			ans.n = 3;
			for (int i = 0; i < 3; i++) {
				ans(i) = getIterator(current->vertex(i));
			}
			
		} else if ( previous->is_in_domain() ) {
		// going through the border edge
			std::vector<VertexHandle> cv = commonVertices(current, previous);
			VertexHandle vertex = vertexHandle(it);
			if ( !isBorder(it) /* from inner node */ ||
				 /* or from border node but first going through inner area */
			     !Utils::has(cv, vertex) ) {
			// return border edge or single border vertex
				ans.n = (int)cv.size();
				for (int i = 0; i < ans.n; i++) {
					ans(i) = getIterator(cv[(size_t)i]);
				}
			}
			
		}
		
		return ans;
	}
	/// @}

	/**
	 * Wrapper around CGAL LineFaceCirculator in order to prevent different
	 * yield situations with line_walk (empty circulator when the line is 
	 * tangent to body and body is at the right of the line, etc..)
	 */
	struct LineWalker {
	
		/** Create line walker from point it in direction determined by shift */
		LineWalker(const Cgal2DGrid* grid_, const Iterator& it, const Real2& shift) :
				grid(grid_), p(grid->coordsD(it)), q(p + shift) {
			assert_true(p != q);

			alongBorder = false;
			if (grid->isBorder(it)) {
				auto borderNeighbors = grid->findBorderNeighbors(it);
				alongBorder = 
						fabs(linal::orientedArea(p, q, grid->coordsD(
								borderNeighbors.first))) < EQUALITY_TOLERANCE ||
						fabs(linal::orientedArea(p, q, grid->coordsD(
								borderNeighbors.second))) < EQUALITY_TOLERANCE;	
			}
			
			// find not empty LineFaceCirculator		
			VertexHandle startVertex = grid->vertexHandle(it);
			plusPlus = true;
			lfc = grid->triangulation.line_walk(cgalPoint2(p), cgalPoint2(q),
					startVertex->incident_faces());
			currentFace = lfc;
			
			if (currentFace == NULL) { // lfc is empty
			// LineFaceCirculator recognizes tangent faces if they are at the left 
			// of the line only; try to change line direction
				plusPlus = false;
				lfc = grid->triangulation.line_walk(cgalPoint2(q), cgalPoint2(p));
				currentFace = lfc;
			}
			
			if (currentFace != NULL) {
			// Now, lfc is not empty.
			// Make it on "in_domain" face which has startVertex
				correctBorderYieldCase();
				auto lfcBegin = lfc;
				while ( !(faceHandle()->is_in_domain() &&
				          faceHandle()->has_vertex(startVertex)) ) {
					next();
					if (lfc == lfcBegin) {
						currentFace = NULL;
						break;
					}
				}
			}
		}
		
		/** Go to the next face and return it */
		FaceHandle next() {
			if (plusPlus) { ++lfc; }
			else          { --lfc; }
			currentFace = lfc;
			correctBorderYieldCase();
			return faceHandle();
		}
		
		FaceHandle faceHandle() const {
			return currentFace;
		}
		bool isValid() const { 
			return currentFace != NULL && currentFace->is_in_domain();
		}
		
	private:
		const Cgal2DGrid * const grid;
		const Real2 p, q; // start and finish
		LineFaceCirculator lfc;
		FaceHandle currentFace; // not always equal to lfc
		
		bool plusPlus; // use "++" not "--" to go to the next face
		bool alongBorder; // lfc goes along the border not inside the body
		
		void correctBorderYieldCase() {
			if (alongBorder && !currentFace->is_in_domain()) {
			// try neighbor face on the other side of the line
				Real2 currentCenter = linal::center({
						real2(currentFace->vertex(0)->point()), 
						real2(currentFace->vertex(1)->point()), 
						real2(currentFace->vertex(2)->point()) });
				for (int i = 0; i < 3; i++) {
					FaceHandle neighbor = currentFace->neighbor(i);
					if ( !neighbor->is_in_domain() ) { continue; }
					Real2 neighborCenter = linal::center({
							real2(neighbor->vertex(0)->point()), 
							real2(neighbor->vertex(1)->point()), 
							real2(neighbor->vertex(2)->point()) });
					
					if ( linal::crossProduct(q - p, currentCenter - p) *
					     linal::crossProduct(q - p, neighborCenter - p) < 0 ) {
					// on the different sides of the line
						auto common = grid->commonVertices(currentFace, neighbor);
						if (fabs(linal::orientedArea(
									p, q, real2(common.at(0)->point()))) < EQUALITY_TOLERANCE &&
						    fabs(linal::orientedArea(
						    		p, q, real2(common.at(1)->point()))) < EQUALITY_TOLERANCE) {
						// their common vertices lie on the line
							currentFace = neighbor;
							break;
						}
					}
				}
			}
		}
		
	};

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
