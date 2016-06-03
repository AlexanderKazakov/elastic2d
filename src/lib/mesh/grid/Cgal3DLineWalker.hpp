#ifndef LIBGCM_CGAL3DLINEWALKER_HPP
#define LIBGCM_CGAL3DLINEWALKER_HPP

#include <lib/mesh/grid/Cgal3DGrid.hpp>


namespace gcm {

/** Struct to go along the line cell-by-cell through Cgal3DGrid */
class Cgal3DLineWalker {
	typedef Cgal3DGrid::Iterator      Iterator;
	typedef Cgal3DGrid::VertexHandle  VertexHandle;
	typedef Cgal3DGrid::CellHandle    CellHandle;
	typedef Cgal3DGrid::CgalPoint3    CgalPoint3;
	typedef Cgal3DGrid::Triangulation Triangulation;
	typedef Triangulation::Edge       Edge;

public:
	/** Create line walker from point it to point it + shift */
	Cgal3DLineWalker(const Cgal3DGrid* const grid_, const Iterator it, const Real3 shift) :
			grid(grid_), 
			q(grid->vertexHandle(it)), 
			p(q->point() + Cgal3DGrid::cgalVector3(shift)) {

		t = findCrossedTetrahedron();
		if (t == NULL) { return; }
		
		u = otherVertex(t, q, q, q);
		v = otherVertex(t, q, q, u);
		w = otherVertex(t, q, u, v);
		if (orientation(u, v, w, q->point()) < 0) { std::swap(u, v); }
	}
	
	/** Current cell */
	CellHandle cell() const {
		return t;
	}
	
	/** Go to the next cell along the line and return it */
	CellHandle next() {
		assert_true(t != NULL);
		assert_false(grid->triangulation.is_infinite(t));
		
		assert_ge(orientation(u, w, v, p), 0); // p is not behind already
		
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
		
		return cell();
	}

private:
	const Cgal3DGrid* const grid;
	const VertexHandle q; ///< start point
	const CgalPoint3 p;   ///< finish point
	CellHandle t;         ///< currentCell
	VertexHandle u, v, w; ///< vertices of the face of current cell 
			///< the line goes through
		
	CellHandle findCrossedTetrahedron() const {
	/// return incident to q finite cell which is crossed by ray (q, p)
	/// or NULL if there isn't such (i.e. crossed cell is infinite)
		
		std::list<CellHandle> fic;
		grid->triangulation.finite_incident_cells(q, std::back_inserter(fic));
		
		for (CellHandle candidate : fic) {
			if ( !grid->isInDomain(candidate) ) { continue; }
			VertexHandle a = otherVertex(candidate, q, q, q);
			VertexHandle b = otherVertex(candidate, q, q, a);
			VertexHandle c = otherVertex(candidate, q, a, b);
			
			Real4 lambda = linal::barycentricCoordinates(
					Cgal3DGrid::real3(q->point()), Cgal3DGrid::real3(a->point()),
					Cgal3DGrid::real3(b->point()), Cgal3DGrid::real3(c->point()),
					Cgal3DGrid::real3(p));
			
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
		for (int i = 0; i < 4; i++) {
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
			const VertexHandle c, const CgalPoint3 d) {
		return linal::orientedVolume(
				Cgal3DGrid::real3(a->point()), Cgal3DGrid::real3(b->point()),
				Cgal3DGrid::real3(c->point()), Cgal3DGrid::real3(d));
	}
	
};


}

#endif // LIBGCM_CGAL3DLINEWALKER_HPP
