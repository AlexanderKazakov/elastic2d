#ifndef LIBGCM_CGALTRIANGULATION_HPP
#define LIBGCM_CGALTRIANGULATION_HPP

#include <libgcm/engine/AbstractEngine.hpp>
#include <libgcm/grid/simplex/cgal/Cgal3DTriangulation.hpp>
#include <libgcm/grid/simplex/cgal/Cgal2DTriangulation.hpp>

namespace gcm {

template<typename, int> class LineWalker;


template<int Dimensionality, typename VertexInfo, typename CellInfo>
struct CgalTriangulationBase;

template<typename VertexInfo, typename CellInfo>
struct CgalTriangulationBase<2, VertexInfo, CellInfo> {
	typedef Cgal2DTriangulation<VertexInfo, CellInfo> type;
};

template<typename VertexInfo, typename CellInfo>
struct CgalTriangulationBase<3, VertexInfo, CellInfo> {
	typedef Cgal3DTriangulation<VertexInfo, CellInfo> type;
};


/**
 * Triangulation in Dimensionality space by CGAL library.
 * @tparam Dimensionality space dimensionality
 * @tparam VertexInfo type of auxiliary information stored in vertices
 * @tparam CellInfo   type of auxiliary information stored in cells
 */
template<int Dimensionality, typename VertexInfo, typename CellInfo>
class CgalTriangulation :
		public CgalTriangulationBase<Dimensionality, VertexInfo, CellInfo>::type {
public:
	
	typedef typename CgalTriangulationBase<
			Dimensionality, VertexInfo, CellInfo>::type Base;
	
	/// Some sort of pointer to triangulation vertex
	typedef typename Base::VertexHandle                 VertexHandle;
	/// Some sort of pointer to triangulation cell
	typedef typename Base::CellHandle                   CellHandle;
	
	/// Point (Vector) in Dimensionality space
	typedef typename Base::RealD                        RealD;
	
	
	/// Space dimensionality
	static const int DIMENSIONALITY = Dimensionality;
	
	/// Number of vertices in cell
	static const int CELL_POINTS_NUMBER = DIMENSIONALITY + 1;
	/// Number of vertices in face
	static const int FACE_POINTS_NUMBER = DIMENSIONALITY;
	
	
	/// An *estimation* of maximal possible number of vertices connected 
	/// with some inner vertex in the grid (it can be more in a very rare cases)
	static const int MAX_NUMBER_OF_NEIGHBOR_VERTICES =
			Base::MAX_NUMBER_OF_NEIGHBOR_VERTICES;
	
	
	/// @name iteration over all finite vertices @{
	typedef typename Base::Triangulation::Finite_vertices_iterator VerticesIterator;
	VerticesIterator verticesBegin() const {
		return this->triangulation.finite_vertices_begin();
	}
	VerticesIterator verticesEnd() const {
		return this->triangulation.finite_vertices_end();
	}
	/// @}
	
	
	CgalTriangulation(const Task& task) : Base(task) {
		rescale(task.simplexGrid.scale);
	}
	virtual ~CgalTriangulation() { }
	
	
	/** Read-only access to points coordinates */
	static RealD coordsD(const VertexHandle vh) {
		return Base::realD(vh->point());
	}
	
	
	/**
	 * Returns common vertices of two given cells.
	 * Optionally, fill in aHasOnly with vertices only the first one has.
	 */
	static std::vector<VertexHandle>
	commonVertices(const CellHandle& a, const CellHandle& b, 
			std::vector<VertexHandle>* aHasOnly = nullptr) {
		
		std::vector<VertexHandle> common;
		for (int i = 0; i < CELL_POINTS_NUMBER; i++) {
			VertexHandle v = a->vertex(i);
			if (b->has_vertex(v)) {
				common.push_back(v);
			} else {
				if (aHasOnly != nullptr) {
					aHasOnly->push_back(v);
				}
			}
		}
		
		return common;
	}
	
	
	/**
	 * Locate cell that contains point on specified distance (shift)
	 * from specified vertex vh by triangulation.locate function.
	 * @threadsafe
	 */
	CellHandle locateOwnerCell(
			const VertexHandle beginVertex, const RealD shift) const {
		auto query = beginVertex->point() + Base::cgalVectorD(shift);
		CellHandle ans;
		#pragma omp critical
		{
		ans = this->triangulation.locate(
				query, Base::someCellOfVertex(beginVertex));
		}
		return ans;
	}
	
	
	/** 
	 * Move specified point on specified distance
	 * without any triangulation reconstruction
	 */
	void move(const VertexHandle vh, const RealD distance) {
		auto& point = vh->point();
		point = point + cgalVectorD(distance);
	}
	
	
	/** Debugging helper */
	void printCell(const CellHandle f, const std::string name) const {
		SUPPRESS_WUNUSED(name);
		LOG_INFO("Cell " << name << ":");
		for (int i = 0; i < CELL_POINTS_NUMBER; i++) {
			if (this->triangulation.is_infinite(f->vertex(i))) {
				LOG_INFO("\nINFINITE\n");
			} else {
				LOG_INFO(Base::realD(f->vertex(i)->point()));
			}
		}
	}
	
	
private:
	USE_AND_INIT_LOGGER("gcm.CgalTriangulation")
	friend class LineWalker<CgalTriangulation, DIMENSIONALITY>;
	
	
	/** Scale the triangulation in space */
	void rescale(const real scale) {
		if (scale == 1) { return; }
		for (auto v  = this->triangulation.finite_vertices_begin();
		          v != this->triangulation.finite_vertices_end(); ++v) {
			auto& point = v->point();
			point = Base::cgalPointD(Base::realD(point) / scale);
		}
	}
	
	
};


}


#endif // LIBGCM_CGALTRIANGULATION_HPP
