#ifndef LIBGCM_CGAL3DGRID_HPP
#define LIBGCM_CGAL3DGRID_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <lib/mesh/grid/UnstructuredGrid.hpp>


namespace gcm {
class Cgal3DMesher;

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
	static const int MAX_NUMBER_OF_NEIGHBOR_VERTICES = 20; // FIXME
	
	
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
	 * from specified point (it) by triangulation.locate function. 
	 * If found face is not "in_domain", returned tetrahedron.valid == false, 
	 * tetrahedron points are not set.
	 * @note for convex bodies result is the same with findOwnerCell
	 */
	Cell locateOwnerCell(const Iterator& it, const Real3& shift) const;
	
	/**
	 * Starting from specified point along the line in specified direction,
	 * find the nearest border face crossed by the line.
	 * @note length of shift must be enough to reach crossing border
	 * @note cases when go from border node outside the body aren't handled
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
	void move(const Iterator& it, const Real3& d) {
		assert_true(movable);
		auto& point = vertexHandle(it)->point();
		point = point + cgalVector3(d);
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
	
	CellHandle findCrossedCell(const VertexHandle start, const Real3& direction) const;
	
	std::vector<VertexHandle> commonVertices(const CellHandle& a, const CellHandle& b,
			std::vector<VertexHandle>* aHasOnly = nullptr) const;

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


	/**
	 * Struct to go along the line cell-by-cell between two points. 
	 */
	struct LineWalker {
		/** Create line walker from coords(it) to coords(it) + shift */
		LineWalker(const Cgal3DGrid* const grid_, const Iterator& it, const Real3& shift_) :
				grid(grid_), start(grid->coords(it)), shift(shift_) {
			
			CellHandle firstCell = grid->findCrossedCell(grid->vertexHandle(it), shift);
			cellsAlongLine.push_back(firstCell);
			if (grid->isInDomain(firstCell)) {
				CellHandle secondCell = grid->triangulation.locate(
						cgalPoint3(start + shift), firstCell);
				if (firstCell != secondCell) {
					cellsAlongLine.push_back(secondCell);
					findCells(0, cellsAlongLine.begin(), 0, 1);
				}
			}
			currentCellIter = cellsAlongLine.begin();
		}
		
		void next() {
			++currentCellIter;
		}
		
		CellHandle cell() const {
			assert_true(currentCellIter != cellsAlongLine.end());
			return correctBorderYieldCase(*currentCellIter);
		}
		
	private:
		const Cgal3DGrid* const grid;
		const Real3 start; ///< start point on the line
		const Real3 shift; ///< vector from start to finish point
		std::list<CellHandle> cellsAlongLine; ///< all cells intersected 
				///< by the segment of the line from start to start + shift
		
		typedef std::list<CellHandle>::iterator CellIterator;
		CellIterator currentCellIter;
		
		void findCells(const int iterationCounter, const CellIterator first, 
				const real firstPosition, const real secondPosition) {
			
			const CellIterator second = std::next(first);
			if ( (*first)->has_neighbor(*second) ) { return; }
			if (iterationCounter > 20) {
				assert_false(grid->isInDomain(*second));
				return; // this is some border case TODO - degenerate case is not handled
			}
			
			real middlePosition = (firstPosition + secondPosition) / 2;
			const CellHandle middleCell = grid->triangulation.locate(
					cgalPoint3(start + shift * middlePosition), *first);
			
			if (middleCell == *first) {
				findCells(iterationCounter + 1, first, middlePosition, secondPosition);
			} else if (middleCell == *second) {
				findCells(iterationCounter + 1, first, firstPosition, middlePosition);
			} else {
				CellIterator middle = cellsAlongLine.insert(second, middleCell);
				findCells(0, first, firstPosition, middlePosition);
				findCells(0, middle, middlePosition, secondPosition);
			}
		}
		
		CellHandle correctBorderYieldCase(const CellHandle c) const {
			return c; // TODO ?
		}
		
		USE_AND_INIT_LOGGER("Cgal3DGridLineWalker")
	};
	

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

	friend class Cgal3DMesher;
	USE_AND_INIT_LOGGER("gcm.Cgal3DGrid")
};


}

#endif // LIBGCM_CGAL3DGRID_HPP
