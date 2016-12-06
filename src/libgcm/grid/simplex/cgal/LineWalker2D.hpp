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
	
	static int otherVertexIndex(
			const CellHandle cell, const VertexHandle a, const VertexHandle b) {
	/// return index of that vertex of cell, which is not a, b
		for (int i = 0; i < CELL_SIZE; i++) {
			VertexHandle d = cell->vertex(i);
			if ( (d != a) && (d != b) ) { return i; }
		}
		THROW_BAD_MESH("Cell contains equal vertices");
	}
	
	static CellHandle neighborThrough(
			const CellHandle cell, const VertexHandle a, const VertexHandle b) {
	/// return neighbor cell that shares with given cell vertices a, b
		return cell->neighbor(otherVertexIndex(cell, a, b));
	}
	
	static VertexHandle otherVertex(
			const CellHandle cell, const VertexHandle a, const VertexHandle b) {
	/// return that vertex of cell, which is not a, b
		return cell->vertex(otherVertexIndex(cell, a, b));
	}
	
	/**
	 * Collect cells along the line from q to p.
	 * The most useful values for q1 and p1 is: q1 = q->point(), p1 = p,
	 * but to handle numerical inexactness along edges and facets,
	 * sometimes it is useful to shift slightly q1 from q or/and p1 from p.
	 * The search accepts only "valid" cells in terms of given predicate:
	 * if meet a not "valid" cell, the search is terminated.
	 */
	template<typename Predicate>
	static std::vector<CellHandle> cellsAlongSegment(
			const Triangulation* triangulation, const Predicate isValid,
			const VertexHandle q, const Real2 p, const Real2 q1, const Real2 p1) {
		
		std::vector<CellHandle> ans;
		CellHandle t = findCrossedCell(triangulation, isValid, q, p);
		if (t == NULL) { return ans; }
		
		VertexHandle l = otherVertex(t, q, q);
		VertexHandle r = otherVertex(t, q, l);
		if (orientation(r, l, q1) < 0) { std::swap(l, r); }
		
		ans.push_back(t);
		while (orientation(p1, r, l) < 0) {
			t = neighborThrough(t, r, l);
			ans.push_back(t);
			if (!isValid(t)) { break; }
			
			VertexHandle s = otherVertex(t, r, l);
			if (orientation(s, q1, p1) < 0) {
				r = s;
			} else {
				l = s;
			}
		}
		return ans;
	}
	
	template<typename Predicate>
	static CellHandle findCrossedCell(
			const Triangulation* triangulation, const Predicate isValid,
			const VertexHandle q, const Real2 p) {
	/// return incident to q "valid" cell which is crossed by the ray (q, p)
	/// or NULL if there isn't such (i.e. crossed cell is not "valid")
		std::list<CellHandle> cells = triangulation->allIncidentCells(q);
		for (CellHandle candidate : cells) {
			if (!isValid(candidate)) { continue; }
			
			VertexHandle a = otherVertex(candidate, q, q);
			VertexHandle b = otherVertex(candidate, q, a);
			Real3 lambda = linal::barycentricCoordinates(
					Triangulation::realD(q), Triangulation::realD(a),
					Triangulation::realD(b), Triangulation::realD(p));
			
			if (lambda(0) <= 1 &&
					lambda(1) >= 0 && lambda(2) >= 0) {
				return candidate;
			}
		}
		return NULL;
	}
	
};

}


#endif // LIBGCM_LINEWALKER2D_HPP
