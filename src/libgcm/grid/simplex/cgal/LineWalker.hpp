#ifndef LIBGCM_LINEWALKER_HPP
#define LIBGCM_LINEWALKER_HPP

#include <libgcm/grid/simplex/cgal/LineWalker2D.hpp>


namespace gcm {

/** @see all comments in LineWalker2D */
template<typename Triangulation>
class LineWalker<Triangulation, 3> {
public:
	typedef typename Triangulation::VertexHandle  VertexHandle;
	typedef typename Triangulation::CellHandle    CellHandle;
	static const int CELL_SIZE = Triangulation::CELL_POINTS_NUMBER;
	LineWalker() {}
	
	template<typename TA, typename TB, typename TC, typename TD>
	static real orientation(const TA a, const TB b, const TC c, const TD d) {
		return linal::orientedVolume(
				Triangulation::realD(a), Triangulation::realD(b),
				Triangulation::realD(c), Triangulation::realD(d));
	}
	
	template<typename Predicate>
	static std::vector<CellHandle> collectCells(
			const Predicate isValid, const Real3 q, const Real3 p,
			CellHandle t, VertexHandle u, VertexHandle v, VertexHandle w) {
		std::vector<CellHandle> ans;
		ans.push_back(t);
		while (orientation(u, v, w, p) < 0) {
			t = Triangulation::neighborThrough(t, u, v, w);
			ans.push_back(t);
			if (!isValid(t)) { break; }
			
			VertexHandle s = Triangulation::otherVertex(t, u, v, w);
			if (orientation(u, s, q, p) > 0) {
				if (orientation(v, s, q, p) > 0) {
					u = s;
				} else {
					w = s;
				}
			} else {
				if (orientation(w, s, q, p) > 0) {
					v = s;
				} else {
					u = s;
				}
			}
		}
		return ans;
	}
	
	template<typename Predicate>
	static std::vector<CellHandle> cellsAlongSegment(
			const Triangulation* triangulation, const Predicate isValid,
			const VertexHandle q, const Real3 p) {
		static_assert(Triangulation::DIMENSIONALITY == 3, "");
		
		CellHandle t = triangulation->findCrossedIncidentCell(isValid, q, p, 0);
		if (t == NULL) { return std::vector<CellHandle>(); }
		
		VertexHandle u = Triangulation::otherVertex(t, q, q, q);
		VertexHandle v = Triangulation::otherVertex(t, q, q, u);
		VertexHandle w = Triangulation::otherVertex(t, q, u, v);
		if (orientation(u, v, w, q) < 0) { std::swap(u, v); }
		
		return collectCells(isValid, Triangulation::realD(q), p, t, u, v, w);
	}
	
	template<typename Predicate>
	static std::vector<CellHandle> cellsAlongSegment(
			const Triangulation* triangulation, const Predicate isValid,
			const CellHandle t, const Real3 p) {
		const Real3 q = Triangulation::center(t);
		
		VertexHandle u = NULL, v = NULL, w = NULL;
		triangulation->findCrossedInsideOutFacet(t, q, p, u, v, w, 0);
		if (u == NULL) { 
			triangulation->findCrossedInsideOutFacet(
					t, q, p, u, v, w, EQUALITY_TOLERANCE);
		}
		if (u == NULL) {
			return std::vector<CellHandle>();
		}
		
		if (orientation(u, v, w, q) < 0) { std::swap(u, v); }
		return collectCells(isValid, q, p, t, u, v, w);
	}
	
};

}

#endif // LIBGCM_LINEWALKER_HPP
