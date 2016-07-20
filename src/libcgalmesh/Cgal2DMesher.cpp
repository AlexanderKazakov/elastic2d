#include <libcgalmesh/Cgal2DMesher.hpp>

#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>


using namespace cgalmesh;


void Cgal2DMesher::
triangulate(const double spatialStep, const std::vector<Body> bodies,
		ResultingTriangulation& result) {
	
	typedef CGAL::Delaunay_mesh_face_base_2<K>                  Fb;
	typedef CGAL::Triangulation_data_structure_2<Vb, Fb>        Tds;
	typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>  CDT;
	typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>            Criteria;
	typedef CGAL::Delaunay_mesher_2<CDT, Criteria>              Mesher;
	
	CDT cdt;
	
	auto insertPolygon = [&](const Polygon& polygon) {
		// insert points and constraints - lines between them - to cdt
		auto point = polygon.vertices_begin();
		CDT::Vertex_handle first = cdt.insert(*point);
		CDT::Vertex_handle last = first;
		++point;
		
		while(point != polygon.vertices_end()) {
			CDT::Vertex_handle current = cdt.insert(*point);
			cdt.insert_constraint(last, current);
			last = current;
			++point;
		}
		
		cdt.insert_constraint(last, first);
	};
	
	// Special "seeds" in the interior of inner cavities
	// to tell CGAL do not mesh these cavities
	std::list<CgalPoint2> listOfSeeds;
	
	// insert all task outer borders and inner cavities
	for (const auto& body : bodies) {
		insertPolygon(makePolygon(body.outer));
		for (const auto& innerCavity : body.inner) {
			insertPolygon(makePolygon(innerCavity));
			listOfSeeds.push_back(findInnerPoint(makePolygon(innerCavity)));
		}
	}
	
	Mesher mesher(cdt);
	mesher.set_seeds(listOfSeeds.begin(), listOfSeeds.end());
	
	Criteria meshingCriteria;
	assert(spatialStep > 0);
	meshingCriteria.set_size_bound(spatialStep);
	mesher.set_criteria(meshingCriteria);
	
	mesher.refine_mesh(); //< meshing
	
	copyTriangulation(cdt, result);
}


Cgal2DMesher::Polygon Cgal2DMesher::
makePolygon(const std::vector<Body::Point>& points) {
	/// convert point set to simple CGAL polygon
	assert(points.size() >= 3); // polygon is a closed line
	
	Polygon polygon;
	for (const auto& p : points) {
		polygon.push_back(CgalPoint2(p[0], p[1]));
	}
	assert(polygon.is_simple()); // has not intersections
	
	return polygon;
}


Cgal2DMesher::CgalPoint2 Cgal2DMesher::
findInnerPoint(const Polygon& polygon) {
	/// find inner point in polygon
	CgalVector2 alongBorder(polygon.vertex(0), polygon.vertex(1));
	CgalVector2 crossBorder = alongBorder.perpendicular(CGAL::Orientation::CLOCKWISE);
	CgalPoint2 middle = CGAL::midpoint(polygon.vertex(0), polygon.vertex(1));
	
	CgalPoint2 innerPoint = middle + crossBorder;
	int n = 1; // iteration number
	while ( !polygon.has_on_bounded_side(innerPoint) ) {
	// watching on the both sides of the edge, getting closer on each iteration
		innerPoint = middle + crossBorder / pow(-2, n);
		assert(n < 20);
		n++;
	}
	
	return innerPoint;
}

