#ifndef LIBGCM_LINEWALKER2D_HPP
#define LIBGCM_LINEWALKER2D_HPP

#include <libgcm/linal/linal.hpp>

namespace gcm {

/**
 * Struct to go along the line cell-by-cell between two points 
 * through a triangulation in order to find cell that contains second point.
 * Implements "line walk in triangulation" strategy.
 * For some degenerate cases points shifting can be implemented
 * to avoid line misses when go exactly along edges and facets.
 * Partially specialized for 2D and 3D cases.
 */
template<typename Triangulation, int Dimensionality>
struct LineWalker;


template<typename Triangulation>
class LineWalker<Triangulation, 2> {
public:
	typedef typename Triangulation::VertexHandle  VertexHandle;
	typedef typename Triangulation::CellHandle    CellHandle;
	static const int CELL_SIZE = Triangulation::CELL_POINTS_NUMBER;
	
	template<typename TA, typename TB, typename TD>
	static real orientation(const TA a, const TB b, const TD d) {
	/// oriented volume to determine where to go next
		return linal::orientedArea(Triangulation::realD(a),
				Triangulation::realD(b), Triangulation::realD(d));
	}
	
	template<typename Predicate>
	static std::vector<CellHandle> collectCells(
			const Predicate isValid, const Real2 q, const Real2 p,
			CellHandle t, VertexHandle r, VertexHandle l) {
	/// given with start parameters, collect "valid" cells along the segment qp;
	/// if found not "valid" cell, end the search
		std::vector<CellHandle> ans;
		ans.push_back(t);
		while (orientation(p, r, l) < 0) {
			t = Triangulation::neighborThrough(t, r, l);
			ans.push_back(t);
			if (!isValid(t)) { break; }
			
			VertexHandle s = Triangulation::otherVertex(t, r, l);
			if (orientation(s, q, p) < 0) {
				r = s;
			} else {
				l = s;
			}
		}
		return ans;
	}
	
	/**
	 * Collect cells along the line from the vertex q to point p.
	 * The search accepts only "valid" cells in terms of given predicate:
	 * if meet a not "valid" cell, the search is stopped.
	 * This method is unstable to numerical inexactness, 
	 * like going along borders and very long lines qp
	 */
	template<typename Predicate>
	static std::vector<CellHandle> cellsAlongSegment(
			const Triangulation* triangulation, const Predicate isValid,
			const VertexHandle q, const Real2 p) {
		static_assert(Triangulation::DIMENSIONALITY == 2, "");
		
		CellHandle t = triangulation->findCrossedIncidentCell(isValid, q, p, 0);
		if (t == NULL) { return std::vector<CellHandle>(); }
		
		VertexHandle l = Triangulation::otherVertex(t, q, q);
		VertexHandle r = Triangulation::otherVertex(t, q, l);
		if (orientation(r, l, q) < 0) { std::swap(l, r); }
		
		return collectCells(isValid, Triangulation::realD(q), p, t, r, l);
	}
	
	/**
	 * Collect cells along the line from q1 == (q + delta) to p.
	 * This is a generalization of the usual algo for the case delta = 0.
	 * In order to handle numerical inexactness along edges and facets,
	 * sometimes it is useful to *slightly* shift q1 from q.
	 * The search accepts only "valid" cells in terms of given predicate:
	 * if meet a not "valid" cell, the search is stopped.
	 */
	template<typename Predicate>
	static std::vector<CellHandle> cellsAlongSegment(
			const Triangulation* triangulation, const Predicate isValid,
			const VertexHandle q, const Real2 p, const Real2 q1) {
		
		CellHandle t = triangulation->findContainingIncidentCell(isValid, q, q1, 0);
		if (t == NULL) { return std::vector<CellHandle>(); }
		
		VertexHandle l = NULL, r = NULL;
		triangulation->findCrossedInsideOutFacet(t, q1, p, l, r, 0);
		if (l == NULL) { return std::vector<CellHandle>(); }
		if (orientation(r, l, q1) < 0) { std::swap(l, r); }
		
		return collectCells(isValid, q1, p, t, r, l);
	}
	
};

}


#endif // LIBGCM_LINEWALKER2D_HPP
