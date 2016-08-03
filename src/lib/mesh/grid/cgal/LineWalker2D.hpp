#ifndef LIBGCM_LINEWALKER2D_HPP
#define LIBGCM_LINEWALKER2D_HPP

#include <lib/linal/linal.hpp>

namespace gcm {

template<typename Triangulation, int Dimensionality>
struct LineWalker;


/** 
 * Struct to go along the line cell-by-cell through a 2D triangulation.
 * Implement "line walk in triangulation" strategy.
 * For some cases like queries along borders
 * infinite loops and line missings are possible.
 */
template<typename Triangulation>
class LineWalker<Triangulation, 2> {
public:
	
	typedef typename Triangulation::VertexHandle  VertexHandle;
	typedef typename Triangulation::CellHandle    CellHandle;
	typedef typename Triangulation::CgalPointD    CgalPointD;
	
	static const int TETR_SIZE = Triangulation::CELL_POINTS_NUMBER;
	
	/**
	 * Create line walker from point vh to point vh + shift.
	 * Line walker recognizes cells with id == id_ only
	 */
	LineWalker(const Triangulation* const triangulation_,
			const VertexHandle vh, const Real2 shift, size_t id_) :
					triangulation(triangulation_), 
					q(vh), 
					p(q->point() + Triangulation::cgalVectorD(shift)),
					id(id_) {
		
		t = findCrossedTetrahedron();
		if (t == NULL) { return; }
		
		l = otherVertex(t, q, q);
		r = otherVertex(t, q, l);
		if (orientation(r, l, q->point()) < 0) { std::swap(l, r); }
	}
	
	
	/** 
	 * Go to the next cell along the line and return it
	 */
	CellHandle next() {
		assert_true(t != NULL);
		assert_false(triangulation->triangulation.is_infinite(t));
		
//		assert_ge(orientation(), 0); // p is not behind already
		
		t = neighborThrough(t, r, l);
		VertexHandle s = otherVertex(t, r, l);
		
		if (orientation(s, q, p) < 0) {
			r = s;
		} else {
			l = s;
		}
		
		return t;
	}
	
	
	CellHandle currentCell() const { return t; }
	
	
private:
	const Triangulation* const triangulation;
	const VertexHandle q; ///< start point
	const CgalPointD p;   ///< finish point
	CellHandle t;         ///< currentCell
	VertexHandle l, r; ///< vertices of the face of current cell 
			///< the line goes through
	const size_t id;
	
	CellHandle findCrossedTetrahedron() const {
	/// return incident to q finite cell which is crossed by ray (q, p)
	/// or NULL if there isn't such (i.e. crossed cell is infinite)
		
		std::list<CellHandle> cells = triangulation->allIncidentCells(q);
		
		for (CellHandle candidate : cells) {
			if ( candidate->info().getGridId() != id ) { continue; }
			VertexHandle a = otherVertex(candidate, q, q);
			VertexHandle b = otherVertex(candidate, q, a);
			
			Real3 lambda = linal::barycentricCoordinates(
					Triangulation::realD(q->point()), Triangulation::realD(a->point()),
					Triangulation::realD(b->point()),
					Triangulation::realD(p));
			
			if (lambda(0) > 1) {
			// behind q
				continue;
			}
			if (lambda(1) >= 0 && lambda(2) >= 0) {
			// inner
				return candidate;
			}
		}
		
		return NULL;
	}
	
	static CellHandle neighborThrough(const CellHandle cell,
			const VertexHandle a, const VertexHandle b) {
	/// return neighbor tetrahedron that shares with t vertices a, b
		return cell->neighbor(otherVertexIndex(cell, a, b));
	}
	
	static int otherVertexIndex(const CellHandle cell,
			const VertexHandle a, const VertexHandle b) {
	/// return index of that vertex of cell, which is not a, b
		for (int i = 0; i < TETR_SIZE; i++) {
			VertexHandle d = cell->vertex(i);
			if ( (d != a) && (d != b) ) { return i; }
		}
		THROW_BAD_MESH("Cell seems to be broken");
	}
	
	static VertexHandle otherVertex(const CellHandle cell,
			const VertexHandle a, const VertexHandle b) {
	/// return vertex of cell, which is not a, b
		return cell->vertex(otherVertexIndex(cell, a, b));
	}
	
	static real orientation(const VertexHandle a, const VertexHandle b,
			const CgalPointD d) {
		return linal::orientedArea(
				Triangulation::realD(a->point()), Triangulation::realD(b->point()),
				Triangulation::realD(d));
	}
	
};



/**
 * Struct to go along the line cell-by-cell through a 2D triangulation.
 * Wrapper around CGAL 2D LineFaceCirculator in order to prevent different
 * yield situations with line_walk (empty circulator when the line is 
 * tangent to body and body is at the right of the line, etc..)
 * FIXME - rewrite with simle algo from the article
 *
template<typename Triangulation>
struct LineWalker<Triangulation, 2> {
	
	typedef typename Triangulation::CellHandle                CellHandle;
	typedef typename Triangulation::VertexHandle              VertexHandle;
	typedef typename Triangulation::Base::LineFaceCirculator  LineFaceCirculator;
	

	LineWalker(const Triangulation * const triangulation_, 
			const VertexHandle vh, const Real2& shift) :
					triangulation(triangulation_),
					p(triangulation->coordsD(vh)),
					q(p + shift) {
		
		assert_true(p != q);
		
		alongBorder = false;
		if (grid->isBorder(it)) {
			std::list<> borderNeighbors = grid->findBorderNeighbors(it);
			alongBorder = 
					fabs(linal::orientedArea(p, q, grid->coordsD(
							borderNeighbors.first))) < EQUALITY_TOLERANCE ||
					fabs(linal::orientedArea(p, q, grid->coordsD(
							borderNeighbors.second))) < EQUALITY_TOLERANCE;
		}
		
		// find not empty LineFaceCirculator
		plusPlus = true;
		lfc = triangulation->triangulation.line_walk(
				Triangulation::cgalPointD(p), Triangulation::cgalPointD(q),
				vh->incident_faces());
		currentFace = lfc;
		
		if (currentFace == NULL) { // lfc is empty
		// LineFaceCirculator recognizes tangent faces if they are at the left 
		// of the line only; try to change line direction
			plusPlus = false;
			lfc = triangulation->triangulation.line_walk(
					Triangulation::cgalPointD(q), Triangulation::cgalPointD(p));
			currentFace = lfc;
		}
		
		if (currentFace != NULL) {
		// Now, lfc is not empty.
		// Make it on "in_domain" face which has vh
			correctBorderYieldCase();
			auto lfcBegin = lfc;
			while ( !(triangulation->isInDomain(CellHandle()) &&
					  CellHandle()->has_vertex(vh)) ) {
				next();
				if (lfc == lfcBegin) {
					currentFace = NULL;
					break;
				}
			}
		}
	}
	
	
	CellHandle next() {
		if (plusPlus) { ++lfc; }
		else          { --lfc; }
		currentFace = lfc;
		correctBorderYieldCase();
		return currentFace;
	}
	
	
	CellHandle currentCell() const { return currentFace; }
	
	
private:
	const Triangulation * const triangulation;
	const Real2 p, q; // start and finish
	LineFaceCirculator lfc;
	CellHandle currentFace; // not always equal to lfc
	
	bool plusPlus; // use "++" not "--" to go to the next face
	bool alongBorder; // lfc goes along the border not inside the body
	
	void correctBorderYieldCase() {
		if ( !alongBorder || grid->isInDomain(currentFace)) { return; }
		
		// try neighbor face on the other side of the line
		Real2 currentCenter = linal::center({
				Triangulation::realD(currentFace->vertex(0)->point()), 
				Triangulation::realD(currentFace->vertex(1)->point()), 
				Triangulation::realD(currentFace->vertex(2)->point()) });
		for (int i = 0; i < 3; i++) {
			CellHandle neighbor = currentFace->neighbor(i);
			if ( !grid->isInDomain(neighbor) ) { continue; }
			Real2 neighborCenter = linal::center({
					Triangulation::realD(neighbor->vertex(0)->point()), 
					Triangulation::realD(neighbor->vertex(1)->point()), 
					Triangulation::realD(neighbor->vertex(2)->point()) });
			
			if ( linal::crossProduct(q - p, currentCenter - p) *
				 linal::crossProduct(q - p, neighborCenter - p) < 0 ) {
			// on the different sides of the line
				auto common = grid->commonVertices(currentFace, neighbor);
				if (fabs(linal::orientedArea(
							p, q, Triangulation::realD(common.at(0)->point()))) <
									EQUALITY_TOLERANCE &&
					fabs(linal::orientedArea(
							p, q, Triangulation::realD(common.at(1)->point()))) <
									EQUALITY_TOLERANCE) {
				// their common vertices lie on the line
					currentFace = neighbor;
					break;
				}
			}
		}
	}
	
};
*/

}


#endif // LIBGCM_LINEWALKER2D_HPP
