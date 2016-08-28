#ifndef LIBGCM_LINEWALKER_HPP
#define LIBGCM_LINEWALKER_HPP

#include <libgcm/grid/simplex/cgal/LineWalker2D.hpp>


namespace gcm {

/** 
 * Struct to go along the line cell-by-cell through a 3D triangulation.
 * Implement "line walk in triangulation" strategy.
 * Yield cases seem to be impossible to handle accurately, so for some cases
 * like queries along borders infinite loops and line missings are possible.
 */
template<typename Triangulation>
class LineWalker<Triangulation, 3> {
public:
	
	typedef typename Triangulation::VertexHandle  VertexHandle;
	typedef typename Triangulation::CellHandle    CellHandle;
	typedef typename Triangulation::CgalPointD    CgalPointD;
	
	static const int TETR_SIZE = Triangulation::CELL_POINTS_NUMBER;
	
	/**
	 * Create line walker from point vh to point vh + shift
	 */
	LineWalker(const Triangulation* const triangulation_,
			const VertexHandle vh, const Real3 shift, const size_t id_) :
					triangulation(triangulation_), 
					q(vh), 
					p(q->point() + Triangulation::cgalVectorD(shift)),
					id(id_) {

		t = findCrossedTetrahedron();
		if (t == NULL) { return; }
		
		u = otherVertex(t, q, q, q);
		v = otherVertex(t, q, q, u);
		w = otherVertex(t, q, u, v);
		if (orientation(u, v, w, q->point()) < 0) { std::swap(u, v); }
	}
	
	/** 
	 * Go to the next cell along the line and return it
	 */
	CellHandle next() {
		assert_true(t != NULL);
		assert_false(triangulation->triangulation.is_infinite(t));
		
//		assert_ge(orientation(u, w, v, p), 0); // p is not behind already
		
		t = neighborThrough(t, u, v, w);
		VertexHandle s = otherVertex(t, u, v, w);
		
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
		
		return t;
	}
	
	
	CellHandle currentCell() const { return t; }
	
	
private:
	const Triangulation* const triangulation;
	const VertexHandle q; ///< start point
	const CgalPointD p;   ///< finish point
	CellHandle t;         ///< currentCell
	VertexHandle u, v, w; ///< vertices of the face of current cell 
			///< the line goes through
	const size_t id;
	
	CellHandle findCrossedTetrahedron() const {
	/// return incident to q finite cell which is crossed by ray (q, p)
	/// or NULL if there isn't such (i.e. crossed cell is infinite)
		
		std::list<CellHandle> cells = triangulation->allIncidentCells(q);
		
		for (CellHandle candidate : cells) {
			if ( candidate->info().getGridId() != id ) { continue; }
			VertexHandle a = otherVertex(candidate, q, q, q);
			VertexHandle b = otherVertex(candidate, q, q, a);
			VertexHandle c = otherVertex(candidate, q, a, b);
			
			Real4 lambda = linal::barycentricCoordinates(
					Triangulation::realD(q->point()), Triangulation::realD(a->point()),
					Triangulation::realD(b->point()), Triangulation::realD(c->point()),
					Triangulation::realD(p));
			
			if (lambda(0) > 1) {
			// behind q
				continue;
			}
			if (lambda(1) >= 0 && lambda(2) >= 0 && lambda(3) >= 0) {
			// inner
				return candidate;
			}
		}
		
		return NULL;
	}
	
	static CellHandle neighborThrough(const CellHandle cell,
			const VertexHandle a, const VertexHandle b, const VertexHandle c) {
	/// return neighbor tetrahedron that shares with t vertices a, b, c
		return cell->neighbor(otherVertexIndex(cell, a, b, c));
	}
	
	static int otherVertexIndex(const CellHandle cell,
			const VertexHandle a, const VertexHandle b, const VertexHandle c) {
	/// return index of that vertex of cell, which is not a, b, c
		for (int i = 0; i < TETR_SIZE; i++) {
			VertexHandle d = cell->vertex(i);
			if ( (d != a) && (d != b) && (d != c) ) { return i; }
		}
		THROW_BAD_MESH("Cell seems to be broken");
	}
	
	static VertexHandle otherVertex(const CellHandle cell,
			const VertexHandle a, const VertexHandle b, const VertexHandle c) {
	/// return vertex of cell, which is not a, b, c
		return cell->vertex(otherVertexIndex(cell, a, b, c));
	}
	
	static real orientation(const VertexHandle a, const VertexHandle b,
			const VertexHandle c, const CgalPointD d) {
		return linal::orientedVolume(
				Triangulation::realD(a->point()), Triangulation::realD(b->point()),
				Triangulation::realD(c->point()), Triangulation::realD(d));
	}
	
};


}

#endif // LIBGCM_LINEWALKER_HPP
