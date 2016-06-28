#ifndef LIBGCM_CGAL3DGRID_HPP
#define LIBGCM_CGAL3DGRID_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <lib/mesh/grid/UnstructuredGrid.hpp>


namespace gcm {
class Cgal3DLineWalker;

/**
 * 3D movable unstructured tetrahedron grid by CGAL library
 */
class Cgal3DGrid : public UnstructuredGrid {
public:
	/// Triangulation types
	typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
	typedef CGAL::Triangulation_vertex_base_with_info_3<size_t, K> Vb;
	typedef CGAL::Triangulation_cell_base_with_info_3<size_t, K>   Cb;
	typedef CGAL::Triangulation_data_structure_3<Vb, Cb>           Tds;
	typedef CGAL::Delaunay_triangulation_3<K, Tds>                 Triangulation;
	
	typedef Triangulation::Vertex                                  VertexT;
	typedef Triangulation::Cell                                    CellT;
	typedef Triangulation::Vertex_handle                           VertexHandle;
	typedef Triangulation::Cell_handle                             CellHandle;
	typedef Triangulation::Geom_traits::Vector_3                   CgalVector3;
	typedef Triangulation::Tetrahedron                             CgalTetrahedron;
	typedef Triangulation::Point                                   CgalPoint3;
	typedef Triangulation::Finite_cells_iterator                   FiniteCellsIterator;

	typedef elements::Tetrahedron<Iterator>                        Cell;
	typedef elements::Triangle<Iterator>                           Face;
	
	/// Space dimensionality
	static const int DIMENSIONALITY = 3;
	/// An *estimation* of maximal possible number of vertices connected 
	/// with some inner vertex (it can be more in a very rare cases)
	static const int MAX_NUMBER_OF_NEIGHBOR_VERTICES = 20;
	
	
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
	 * Iteration over all real mesh cells.
	 * "Real" means cell which covers the body not outer space ("isInDomain").
	 * There are also not meshing ("!isInDomain") cells in the triangulation.
	 * They cover the space around bodies, because CGAL triangulation 
	 * is a convex hull of its vertices.
	 */
	///@{
	struct RealCellTester {
		bool operator()(const FiniteCellsIterator& it) const {
			return !grid->isInDomain(it);
		}
		RealCellTester(const Cgal3DGrid* const grid_) : grid(grid_) { }
	private:
		const Cgal3DGrid* const grid;	
	};
	
	typedef CGAL::Filter_iterator<FiniteCellsIterator, RealCellTester> CellIterator;
	
	CellIterator cellBegin() const { 
		return CellIterator(triangulation.finite_cells_end(),
		                    RealCellTester(this),
		                    triangulation.finite_cells_begin());
	}
	CellIterator cellEnd() const { 
		return CellIterator(triangulation.finite_cells_end(),
		                    RealCellTester(this),
		                    triangulation.finite_cells_end());
	}
	///@}
	///@}

	Cgal3DGrid(const Task& task);
	virtual ~Cgal3DGrid() { }


	/** Read-only access to real coordinates */
	Real3 coords(const Iterator& it) const {
		return coordsD(it);
	}
	
	/** Read-only access to real coordinates */
	Real3 coordsD(const Iterator& it) const {
		return real3(vertexHandle(it)->point());
	}
	
	/** Border normal in specified point */
	Real3 normal(const Iterator& it) const;
	
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
	std::array<size_t, 4> getVerticesOfCell(const CellIterator& it) const {
		std::array<size_t, 4> vertices;
		for (size_t i = 0; i < 4; i++) {
			vertices[i] = it->vertex((int)i)->info();
		}
		return vertices;
	}
	
	Iterator getIterator(const VertexHandle v) const {
		return v->info();
	}
	
	size_t cellInfo(const CellHandle c) const {
		return c->info();
	}
	
	bool isInDomain(const CellHandle c) const {
		return ( !triangulation.is_infinite(c) ) && (cellInfo(c) != 0);
	}
	
	bool isBorder(const Iterator& it) const {
		VertexHandle vh = vertexHandle(it);
		std::list<CellHandle> incidentCells;
		triangulation.incident_cells(vh, std::back_inserter(incidentCells));
		for (const auto cell : incidentCells) {
			if (!isInDomain(cell)) {
				return true;
			}
		}
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
	 * Find tetrahedron that contains point on specified distance (shift) 
	 * from specified point (it) by line walk from it to (it+shift).
	 * In some bad degenerate cases cell.n == 0 and points are not set.
	 * If the line goes out of the body immediately (from border in outer 
	 * direction), cell.n == 0 and points are not set.
	 * If the line goes through the body, but go out before reaching
	 * the target point, returned cell.n == 3, and only 3 points are set
	 * - they are three border vertices - the face crossed by the line.
	 * If the line goes through the body, and reach the target point,
	 * returned cell.n == 4, and all 4 points are set.
	 * @note for convex bodies and (n == 4) result is the same with locateOwnerCell
	 */
	Cell findOwnerCell(const Iterator& it, const Real3& shift) const;
	
	/**
	 * Locate tetrahedron that contains point on specified distance (shift) 
	 * from specified point (it) by triangulation.locate function. 
	 * If found cell is not "in_domain" (out of body), returned cell.n == 0, 
	 * and points are not set. Else n == 4 and all 4 points are set.
	 * @note for convex bodies and (n == 4) result is the same with findOwnerCell
	 */
	Cell locateOwnerCell(const Iterator& it, const Real3& shift) const;
	
	/**
	 * @return set of border vertices incident to given vertex
	 */
	std::set<Iterator> findBorderNeighbors(const Iterator& it) const;
		
	/** 
	 * Find vertex with specified coordinates. Throw Exception if there isn't such.
	 */
	Iterator findVertexByCoordinates(const Real3& coordinates) const;


protected:
	/// Data
	///@{
	Triangulation triangulation;             ///< CGAL triangulation data structure
	std::vector<VertexHandle> vertexHandles; ///< CGAL-"pointers" to each grid vertex
	std::vector<size_t> borderIndices;       ///< indices of border vertices in vertexHandles
	std::vector<size_t> innerIndices;        ///< indices of inner vertices in vertexHandles
	real effectiveSpatialStep = 0;           ///< used in triangulation criteria and Courant condition
	bool movable = false;                    ///< deformable(true) or immutable(false) grid
	///@}

	/** Move specified point on specified distance */
	void move(const Iterator& it, const Real3& distance) {
		assert_true(movable);
		auto& point = vertexHandle(it)->point();
		point = point + cgalVector3(distance);
	}

	
private:
	/** Functions for building the triangulation */
	///@{
	void markInnersAndBorders();
	///@}

	/// Auxilliary functions
	///@{
	
	VertexHandle vertexHandle(const Iterator it) const {
		return vertexHandles[getIndex(it)];
	}
	
	std::vector<VertexHandle> commonVertices(const CellHandle& a, const CellHandle& b,
			std::vector<VertexHandle>* aHasOnly = nullptr) const;
	
	bool contains(const CellHandle cell, const Real3& q) const {
	/// is tetrahedron with small layer around contains the point 
		Real3 a = real3(cell->vertex(0)->point());
		Real3 b = real3(cell->vertex(1)->point());
		Real3 c = real3(cell->vertex(2)->point());
		Real3 d = real3(cell->vertex(3)->point());
		Real4 lambda = linal::barycentricCoordinates(a, b, c, d, q);
		return (lambda(0) > -EQUALITY_TOLERANCE) &&
		       (lambda(1) > -EQUALITY_TOLERANCE) &&
		       (lambda(2) > -EQUALITY_TOLERANCE) &&
		       (lambda(3) > -EQUALITY_TOLERANCE);
	}

	Cell createCell(const Iterator& it,
			const CellHandle current, const CellHandle previous) const {
	/// create Cell used as answer to numerical method queries about point location
		Cell ans;
		ans.n = 0;
		if (current == NULL) {
			assert_true(previous == NULL);
			return ans;
		}
		
		if (isInDomain(current)) {
			assert_true(isInDomain(previous));
			ans.n = 4;
			for (int i = 0; i < 4; i++) {
				ans(i) = getIterator(current->vertex(i));
			}
			
		} else if ( isInDomain(previous) ) {
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
	///@}

	void printCell(const CellHandle& f, const std::string& name = "") const;

	USE_AND_INIT_LOGGER("gcm.Cgal3DGrid")
	
	friend class Cgal3DLineWalker;

	
public:
	/// @name convertion between CGAL and gcm data types @{
	static CgalPoint3 cgalPoint3(const Real3& p) {
		return CgalPoint3(p(0), p(1), p(2));
	}
	
	static Real3 real3(const CgalPoint3& p) {
		return {p.x(), p.y(), p.z()};
	}
	
	static CgalVector3 cgalVector3(const Real3& p) {
		return CgalVector3(p(0), p(1), p(2));
	}
	/// @}
};


}

#endif // LIBGCM_CGAL3DGRID_HPP
