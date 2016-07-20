#ifndef LIBGCM_CGAL2DLINEWALKER_HPP
#define LIBGCM_CGAL2DLINEWALKER_HPP

#include <lib/mesh/grid/cgal/Cgal2DGrid.hpp>


namespace gcm {

/**
 * Struct to go along the line cell-by-cell through Cgal2DGrid.
 * Wrapper around CGAL 2D LineFaceCirculator in order to prevent different
 * yield situations with line_walk (empty circulator when the line is 
 * tangent to body and body is at the right of the line, etc..)
 */
struct Cgal2DLineWalker {
	typedef Cgal2DGrid::Iterator                     Iterator;
	typedef Cgal2DGrid::FaceHandle                   FaceHandle;
	typedef Cgal2DGrid::VertexHandle                 VertexHandle;
	typedef Cgal2DGrid::Triangulation                Tr;
	typedef Tr::Line_face_circulator                 LineFaceCirculator;

	/** Create line walker from point it in direction determined by shift */
	Cgal2DLineWalker(const Cgal2DGrid* grid_, const Iterator& it, const Real2& shift) :
			grid(grid_), p(grid->coordsD(it)), q(p + shift) {
		assert_true(p != q);

		alongBorder = false;
		if (grid->isBorder(it)) {
			auto borderNeighbors = grid->findBorderNeighbors(it);
			alongBorder = 
					fabs(linal::orientedArea(p, q, grid->coordsD(
							borderNeighbors.first))) < EQUALITY_TOLERANCE ||
					fabs(linal::orientedArea(p, q, grid->coordsD(
							borderNeighbors.second))) < EQUALITY_TOLERANCE;	
		}
		
		// find not empty LineFaceCirculator		
		VertexHandle startVertex = grid->vertexHandle(it);
		plusPlus = true;
		lfc = grid->triangulation.line_walk(
				Cgal2DGrid::cgalPoint2(p), Cgal2DGrid::cgalPoint2(q),
				startVertex->incident_faces());
		currentFace = lfc;
		
		if (currentFace == NULL) { // lfc is empty
		// LineFaceCirculator recognizes tangent faces if they are at the left 
		// of the line only; try to change line direction
			plusPlus = false;
			lfc = grid->triangulation.line_walk(
					Cgal2DGrid::cgalPoint2(q), Cgal2DGrid::cgalPoint2(p));
			currentFace = lfc;
		}
		
		if (currentFace != NULL) {
		// Now, lfc is not empty.
		// Make it on "in_domain" face which has startVertex
			correctBorderYieldCase();
			auto lfcBegin = lfc;
			while ( !(grid->isInDomain(faceHandle()) &&
					  faceHandle()->has_vertex(startVertex)) ) {
				next();
				if (lfc == lfcBegin) {
					currentFace = NULL;
					break;
				}
			}
		}
	}
	
	/** Go to the next face and return it */
	FaceHandle next() {
		if (plusPlus) { ++lfc; }
		else          { --lfc; }
		currentFace = lfc;
		correctBorderYieldCase();
		return faceHandle();
	}
	
	FaceHandle faceHandle() const {
		return currentFace;
	}
	bool isValid() const { 
		return currentFace != NULL && grid->isInDomain(currentFace);
	}
	
private:
	const Cgal2DGrid * const grid;
	const Real2 p, q; // start and finish
	LineFaceCirculator lfc;
	FaceHandle currentFace; // not always equal to lfc
	
	bool plusPlus; // use "++" not "--" to go to the next face
	bool alongBorder; // lfc goes along the border not inside the body
	
	void correctBorderYieldCase() {
		if ( !alongBorder || grid->isInDomain(currentFace)) { return; }
		
		// try neighbor face on the other side of the line
		Real2 currentCenter = linal::center({
				Cgal2DGrid::real2(currentFace->vertex(0)->point()), 
				Cgal2DGrid::real2(currentFace->vertex(1)->point()), 
				Cgal2DGrid::real2(currentFace->vertex(2)->point()) });
		for (int i = 0; i < 3; i++) {
			FaceHandle neighbor = currentFace->neighbor(i);
			if ( !grid->isInDomain(neighbor) ) { continue; }
			Real2 neighborCenter = linal::center({
					Cgal2DGrid::real2(neighbor->vertex(0)->point()), 
					Cgal2DGrid::real2(neighbor->vertex(1)->point()), 
					Cgal2DGrid::real2(neighbor->vertex(2)->point()) });
			
			if ( linal::crossProduct(q - p, currentCenter - p) *
				 linal::crossProduct(q - p, neighborCenter - p) < 0 ) {
			// on the different sides of the line
				auto common = grid->commonVertices(currentFace, neighbor);
				if (fabs(linal::orientedArea(
							p, q, Cgal2DGrid::real2(common.at(0)->point()))) <
									EQUALITY_TOLERANCE &&
					fabs(linal::orientedArea(
							p, q, Cgal2DGrid::real2(common.at(1)->point()))) <
									EQUALITY_TOLERANCE) {
				// their common vertices lie on the line
					currentFace = neighbor;
					break;
				}
			}
		}
	}
	
};


}


#endif // LIBGCM_CGAL2DLINEWALKER_HPP
