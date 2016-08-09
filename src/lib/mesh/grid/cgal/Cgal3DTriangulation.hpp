#ifndef LIBGCM_CGAL3DTRIANGULATION_HPP
#define LIBGCM_CGAL3DTRIANGULATION_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

#include <libcgalmesher/Cgal3DMesher.hpp>
#include <lib/mesh/mesh_loaders/InmMeshLoader.hpp>

#include <lib/linal/linal.hpp>
#include <lib/util/task/Task.hpp>
#include <lib/util/Logging.hpp>

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
	 */
	std::list<CellHandle> allIncidentCells(const VertexHandle vh) const {
		std::list<CellHandle> ans;
		triangulation.incident_cells(vh, std::back_inserter(ans));
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
	
	static RealD realD(const CgalPointD& p) {
		return {p.x(), p.y(), p.z()};
	}
	
	static CgalVectorD cgalVectorD(const RealD& p) {
		return CgalVectorD(p(0), p(1), p(2));
	}
	/// @}
	
	
	/** Is tetrahedron with a small layer around contains the point */
	static bool contains(const CellHandle cell, const RealD& q) {
		RealD a = realD(cell->vertex(0)->point());
		RealD b = realD(cell->vertex(1)->point());
		RealD c = realD(cell->vertex(2)->point());
		RealD d = realD(cell->vertex(3)->point());
		auto lambda = linal::barycentricCoordinates(a, b, c, d, q);
		return (lambda(0) > -EQUALITY_TOLERANCE) &&
		       (lambda(1) > -EQUALITY_TOLERANCE) &&
		       (lambda(2) > -EQUALITY_TOLERANCE) &&
		       (lambda(3) > -EQUALITY_TOLERANCE);
	}
	
	
protected:
	Triangulation triangulation; ///< CGAL triangulation structure
	
	USE_AND_INIT_LOGGER("gcm.Cgal3DTriangulation")
	friend class CgalTriangulation<DIMENSIONALITY, VertexInfo, CellInfo>;
};


}

#endif // LIBGCM_CGAL3DTRIANGULATION_HPP
