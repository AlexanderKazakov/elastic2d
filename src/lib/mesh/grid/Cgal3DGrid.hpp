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
	 * Find tetrahedron that contains point on specified distance (shift) 
	 * from specified point (it) by line walk from it to (it+shift).
	 * If the line goes out of the body, returned tetrahedron.valid == false, 
	 * tetrahedron points are not set.
	 * @note for convex bodies result is the same with locateOwnerCell
	 */
	Cell findOwnerCell(const Iterator& it, const Real3& shift) const;
	
	/**
	 * Find tetrahedron that contains point on specified distance (shift) 
	 * from specified point (it) by triangulation.locate function with 
	 * correction of border "out of domain" cases.
	 * If found face is not "in_domain", returned tetrahedron.valid == false, 
	 * tetrahedron points are not set.
	 * @note for convex bodies result is the same with findOwnerCell
	 */
	Cell locateOwnerCell(const Iterator& it, const Real3& shift) const;
	
	/**
	 * Starting from specified point along the line in specified direction,
	 * find the nearest border face crossed by the line.
	 * @note length of shift must be enough to reach crossing border
	 * @return crossed border face
	 */
	Face findCrossingBorder(
			const Iterator& start, const Real3& shift) const;
	
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
	std::set<size_t> borderIndices;          ///< indices of border vertices in vertexHandles
	std::set<size_t> innerIndices;           ///< indices of inner vertices in vertexHandles
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
	
	CellHandle correctBorderYieldCase(const CellHandle c) const;
	
	std::set<CellHandle> cellsAround(const CellHandle cell, const int depth = 1) const;
	
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

	Cell createTetrahedron(const CellHandle ch) const {
		Cell ans;
		ans.valid = false;
		
		if (ch == NULL) { return ans; }
		
		if (isInDomain(ch)) {
			ans.valid = true;
			for (int i = 0; i < 4; i++) {
				ans(i) = getIterator(ch->vertex(i));
			}
		}
		
		return ans;
	}
	///@}

	void printCell(const CellHandle& f, const std::string& name = "") const;
 
	static CgalPoint3 cgalPoint3(const Real3& p) {
		return CgalPoint3(p(0), p(1), p(2));
	}
	
	static Real3 real3(const CgalPoint3& p) {
		return {p.x(), p.y(), p.z()};
	}
	
	static CgalVector3 cgalVector3(const Real3& p) {
		return CgalVector3(p(0), p(1), p(2));
	}

	USE_AND_INIT_LOGGER("gcm.Cgal3DGrid")
	
	friend class Cgal3DLineWalker;
};


}

#endif // LIBGCM_CGAL3DGRID_HPP
