#ifndef LIBGCM_CGAL3DTRIANGULATION_HPP
#define LIBGCM_CGAL3DTRIANGULATION_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

#include <libcgalmesher/Cgal3DMesher.hpp>
#include <libgcm/grid/simplex/mesh_loaders/InmMeshLoader.hpp>

#include <libgcm/linal/linal.hpp>
#include <libgcm/util/task/Task.hpp>
#include <libgcm/util/infrastructure/infrastructure.hpp>

namespace gcm {

template<int Dimensionality, typename VertexInfo, typename CellInfo>
class CgalTriangulation;

/**
 * 3D triangulation by CGAL library.
 * Tetrahedron named "cell".
 * @tparam VertexInfo type of auxiliary information stored in vertices
 * @tparam CellInfo   type of auxiliary information stored in cells
 */
template<typename VertexInfo, typename CellInfo>
class Cgal3DTriangulation {
public:
	
	typedef CGAL::Exact_predicates_inexact_constructions_kernel        K;
	typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K> Vb;
	typedef CGAL::Triangulation_cell_base_with_info_3<CellInfo, K>     Cb;
	typedef CGAL::Triangulation_data_structure_3<Vb, Cb>               Tds;
	/// The triangulation type
	typedef CGAL::Delaunay_triangulation_3<K, Tds>                     Triangulation;
	
	typedef typename Triangulation::Vertex_handle           VertexHandle;
	typedef typename Triangulation::Cell_handle             CellHandle;
	typedef typename Triangulation::Geom_traits::Vector_3   CgalVectorD;
	typedef typename Triangulation::Point                   CgalPointD;
	typedef typename Triangulation::All_cells_iterator      AllCellsIterator;
	
	/// Space dimensionality
	static const int DIMENSIONALITY = 3;
	static const int CELL_SIZE = DIMENSIONALITY + 1;
	/// Number of verices in tetrahedra
	static const int TETR_SIZE = DIMENSIONALITY + 1;
	/// Point (Vector) in DIMENSIONALITY space
	typedef linal::Vector<DIMENSIONALITY> RealD;
	/// An *estimation* of maximal possible number of vertices connected 
	/// with some inner vertex in the grid (it can be more in a very rare cases)
	static const int MAX_NUMBER_OF_NEIGHBOR_VERTICES = 20;
	
	
	Cgal3DTriangulation(const Task& task) {
		switch (task.simplexGrid.mesher) {
			case Task::SimplexGrid::Mesher::CGAL_MESHER:
				LOG_DEBUG("Call Cgal3DMesher");
				cgalmesher::Cgal3DMesher::triangulate(
						task.simplexGrid.spatialStep, task.simplexGrid.detectSharpEdges,
						task.simplexGrid.fileName, triangulation);
				break;
			case Task::SimplexGrid::Mesher::INM_MESHER:
				LOG_DEBUG("Call InmMeshLoader");
				InmMeshLoader::load(task.simplexGrid.fileName, triangulation);
				break;
			default:
				THROW_UNSUPPORTED("Unknown mesher");
		}
		
		LOG_INFO("Number of all vertices after meshing: " << triangulation.number_of_vertices());
		LOG_INFO("Number of all cells after meshing: " << triangulation.number_of_cells());
	}
	
	
	/** All cells iteration begin */
	AllCellsIterator allCellsBegin() const {
		return triangulation.all_cells_begin();
	}
	/** All cells iteration end */
	AllCellsIterator allCellsEnd() const {
		return triangulation.all_cells_end();
	}
	
	
	/**
	 * Returns all incident to vh cells. Order of cells is not defined.
	 * @threadsafe
	 */
	std::list<CellHandle> allIncidentCells(const VertexHandle vh) const {
		
		std::list<CellHandle> ans;
		/// Don't know why and how, but reading incident cells from
		/// CGAL triangulation by several threads is not thread-safe.
		/// So we have critical section here
		#pragma omp critical
		{
		triangulation.incident_cells(vh, std::back_inserter(ans));
		}
		return ans;
	}
	
	
	/**
	 * Returns unity normal to contact surface between given cells
	 * (must be neighbors). Direction of normal is from "from" to "to"
	 */
	static RealD contactNormal(const CellHandle from, const CellHandle to) {
		int oppositeVertexIndex = from->index(to);
		VertexHandle opposite = from->vertex(oppositeVertexIndex);
		VertexHandle a = from->vertex((oppositeVertexIndex + 1) % TETR_SIZE);
		VertexHandle b = from->vertex((oppositeVertexIndex + 2) % TETR_SIZE);
		VertexHandle c = from->vertex((oppositeVertexIndex + 3) % TETR_SIZE);
		
		return linal::oppositeFaceNormal(realD(opposite->point()),
				realD(a->point()), realD(b->point()), realD(c->point()));
	}
	
	
	static CellHandle someCellOfVertex(const VertexHandle vh) {
		return vh->cell();
	}
	
	
	static real minimalCellHeight(const CellHandle ch) {
		return linal::minimalHeight(
				realD(ch->vertex(0)->point()),
				realD(ch->vertex(1)->point()),
				realD(ch->vertex(2)->point()),
				realD(ch->vertex(3)->point()));
	}
	
	
	/// @name convertion between CGAL and gcm data types
	/// @{
	static CgalPointD cgalPointD(const RealD& p) {
		return CgalPointD(p(0), p(1), p(2));
	}
	
	static RealD realD(const VertexHandle vh) {
		return realD(vh->point());
	}
	
	static RealD realD(const CgalPointD& p) {
		return {p.x(), p.y(), p.z()};
	}
	
	static RealD realD(const Real3 r) {
		return r;
	}
	
	static CgalVectorD cgalVectorD(const RealD& p) {
		return CgalVectorD(p(0), p(1), p(2));
	}
	/// @}
	
	
	/** Is the cell with a small layer around contains the point */
	static bool contains(const CellHandle cell, const RealD& q, const real eps) {
		RealD a = realD(cell->vertex(0));
		RealD b = realD(cell->vertex(1));
		RealD c = realD(cell->vertex(2));
		RealD d = realD(cell->vertex(3));
		return linal::tetrahedronContains(a, b, c, d, q, eps);
	}
	
	
	/**
	 * Returns incident to vh "valid" cell which is crossed by the ray 
	 * from vh to query or NULL if there isn't such 
	 * (i.e. crossed cell is not "valid")
	 */
	template<typename Predicate>
	CellHandle findCrossedIncidentCell(const Predicate isValid,
			const VertexHandle vh, const Real3 query, const real eps) const {
		std::list<CellHandle> cells = allIncidentCells(vh);
		for (CellHandle candidate : cells) {
			if (!isValid(candidate)) { continue; }
			
			VertexHandle a = otherVertex(candidate, vh, vh, vh);
			VertexHandle b = otherVertex(candidate, vh, vh, a);
			VertexHandle c = otherVertex(candidate, vh, a, b);
			
			if (linal::solidAngleContains(
					realD(vh), realD(a), realD(b), realD(c), query, eps)) {
				return candidate;
			}
		}
		return NULL;
	}
	
	
	/**
	 * The point q must lie inside the tetrahedron t.
	 * Find the face of t which is crossed by the ray qp.
	 * Write the result as a triple of vertices to a,b,c.
	 */
	static void findCrossedInsideOutFacet(
			const CellHandle t, const Real3 q, const Real3 p,
			VertexHandle& a, VertexHandle& b, VertexHandle& c, const real eps) {
		a = NULL; b = NULL; c = NULL;
		for (int i = 0; i < CELL_SIZE; i++) {
			VertexHandle a1 = t->vertex((i + 1) % CELL_SIZE);
			VertexHandle b1 = t->vertex((i + 2) % CELL_SIZE);
			VertexHandle c1 = t->vertex((i + 3) % CELL_SIZE);
			if (!linal::isDegenerate(q, realD(a1), realD(b1), realD(c1), EQUALITY_TOLERANCE) &&
					linal::solidAngleContains(q, realD(a1), realD(b1), realD(c1), p, eps)) {
				a = a1; b = b1; c = c1; break;
			}
		}
	}
	
	
	static int otherVertexIndex(const CellHandle cell,
			const VertexHandle a, const VertexHandle b, const VertexHandle c) {
	/// return index of that vertex of cell, which is not a, b, c
		for (int i = 0; i < CELL_SIZE; i++) {
			VertexHandle d = cell->vertex(i);
			if ( (d != a) && (d != b) && (d != c) ) { return i; }
		}
		THROW_BAD_MESH("Cell contains equal vertices");
	}
	
	static CellHandle neighborThrough(const CellHandle cell,
			const VertexHandle a, const VertexHandle b, const VertexHandle c) {
	/// return neighbor cell that shares with given cell vertices a, b, c
		return cell->neighbor(otherVertexIndex(cell, a, b, c));
	}
	
	static VertexHandle otherVertex(const CellHandle cell,
			const VertexHandle a, const VertexHandle b, const VertexHandle c) {
	/// return that vertex of cell, which is not a, b, c
		return cell->vertex(otherVertexIndex(cell, a, b, c));
	}
	
	
protected:
	Triangulation triangulation; ///< CGAL triangulation structure
	
	USE_AND_INIT_LOGGER("gcm.Cgal3DTriangulation")
	friend class CgalTriangulation<DIMENSIONALITY, VertexInfo, CellInfo>;
};


}

#endif // LIBGCM_CGAL3DTRIANGULATION_HPP
