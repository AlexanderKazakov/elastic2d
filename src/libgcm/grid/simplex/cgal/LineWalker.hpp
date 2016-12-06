#ifndef LIBGCM_LINEWALKER_HPP
#define LIBGCM_LINEWALKER_HPP

#include <libgcm/grid/simplex/cgal/LineWalker2D.hpp>


namespace gcm {

template<typename Triangulation>
class LineWalker<Triangulation, 3> {
public:
	typedef typename Triangulation::VertexHandle  VertexHandle;
	typedef typename Triangulation::CellHandle    CellHandle;
	static const int CELL_SIZE = Triangulation::CELL_POINTS_NUMBER;
	
	template<typename TA, typename TB, typename TC, typename TD>
	static real orientation(const TA a, const TB b, const TC c, const TD d) {
	/// oriented volume to determine where to go next
		return linal::orientedVolume(
				Triangulation::realD(a), Triangulation::realD(b),
				Triangulation::realD(c), Triangulation::realD(d));
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
	
	/** @see comments in LineWalker2D */
	template<typename Predicate>
	static std::vector<CellHandle> cellsAlongSegment(
			const Triangulation* triangulation, const Predicate isValid,
			const VertexHandle q, const Real3 p, const Real3 q1, const Real3 p1) {
		
		std::vector<CellHandle> ans;
		CellHandle t = findCrossedCell(triangulation, isValid, q, p);
		if (t == NULL) { return ans; }
		
		VertexHandle u = otherVertex(t, q, q, q);
		VertexHandle v = otherVertex(t, q, q, u);
		VertexHandle w = otherVertex(t, q, u, v);
		if (orientation(u, v, w, q1) < 0) { std::swap(u, v); }
		
		ans.push_back(t);
		while (orientation(u, v, w, p1) < 0) {
			t = neighborThrough(t, u, v, w);
			ans.push_back(t);
			if (!isValid(t)) { break; }
			
			VertexHandle s = otherVertex(t, u, v, w);
			if (orientation(u, s, q1, p1) > 0) {
				if (orientation(v, s, q1, p1) > 0) {
					u = s;
				} else {
					w = s;
				}
			} else {
				if (orientation(w, s, q1, p1) > 0) {
					v = s;
				} else {
					u = s;
				}
			}
		}
		return ans;
	}
	
	template<typename Predicate>
	static CellHandle findCrossedCell(
			const Triangulation* triangulation, const Predicate isValid,
			const VertexHandle q, const Real3 p) {
		/// return incident to q "valid" cell which is crossed by the ray (q, p)
		/// or NULL if there isn't such (i.e. crossed cell is not "valid")
		std::list<CellHandle> cells = triangulation->allIncidentCells(q);
		for (CellHandle candidate : cells) {
			if (!isValid(candidate)) { continue; }
			
			VertexHandle a = otherVertex(candidate, q, q, q);
			VertexHandle b = otherVertex(candidate, q, q, a);
			VertexHandle c = otherVertex(candidate, q, a, b);
			Real4 lambda = linal::barycentricCoordinates(
					Triangulation::realD(q), Triangulation::realD(a),
					Triangulation::realD(b), Triangulation::realD(c),
					Triangulation::realD(p));
			
			if (lambda(0) <= 1 &&
					lambda(1) >= 0 && lambda(2) >= 0 && lambda(3) >= 0) {
				return candidate;
			}
		}
		return NULL;
	}
	
};

}

#endif // LIBGCM_LINEWALKER_HPP
