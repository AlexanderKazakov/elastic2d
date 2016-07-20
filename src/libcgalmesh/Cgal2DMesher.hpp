#ifndef LIBCGALMESH_CGAL2DMESHER_HPP
#define LIBCGALMESH_CGAL2DMESHER_HPP

#include <string>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Polygon_2.h>


namespace cgalmesh {

/**
 * 2D mesher by CGAL library
 */
class Cgal2DMesher {
public:
	typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
	typedef CGAL::Triangulation_vertex_base_with_info_2<size_t, K> Vb;
	typedef CGAL::Triangulation_face_base_with_info_2<size_t, K>   Cb;
	typedef CGAL::Triangulation_data_structure_2<Vb, Cb>           Tds;
	typedef CGAL::Delaunay_triangulation_2<K, Tds>                 ResultingTriangulation;
	
	typedef ResultingTriangulation::Point                          CgalPoint2;
	typedef ResultingTriangulation::Geom_traits::Vector_2          CgalVector2;
	typedef CGAL::Polygon_2<K, std::vector<CgalPoint2>>            Polygon;
	
	typedef ResultingTriangulation::Vertex                         ResultingVertex;
	typedef ResultingTriangulation::Face                           ResultingCell;
	
	struct Body {
	/// closed polygon without self-intersections
		typedef std::array<double, 2> Point;
		typedef std::vector<Point>    Border;
		
		Border outer;               ///< outer border of the body
		std::vector<Border> inner;  ///< borders of the inner cavities
	};
	
	/**
	 * Build the grid on given geometry
	 * @param spatialStep         effective spatial step
	 * @param bodies              list of bodies to construct
	 * @param result              triangulation to write the result in
	 */
	static void triangulate(
			const double spatialStep, const std::vector<Body> bodies,
			ResultingTriangulation& result);
	
	
private:
	static Polygon makePolygon(const Body::Border& points);
	static CgalPoint2 findInnerPoint(const Polygon& polygon);
	
	
	/** Load result triangulation from mesher triangulation */
	template<typename MeshingTriangulation>
	static void copyTriangulation(const MeshingTriangulation& mt, 
			ResultingTriangulation& result) {
		
		typedef typename MeshingTriangulation::Vertex    MeshingVertex;
		typedef typename MeshingTriangulation::Face      MeshingCell;
		
		// converter structures for CGAL copy_tds function
		struct VertexConverter {
			ResultingVertex operator()(const MeshingVertex& src) const {
				return ResultingVertex(src.point());
			}
			void operator()(const MeshingVertex&, ResultingVertex&) const { }
		} vertexConverter;
		struct CellConverter {
			ResultingCell operator()(const MeshingCell& src) const {
				ResultingCell cellT;
				cellT.info() = src.is_in_domain() ? 1 : 0;
				return cellT;
			}
			void operator()(const MeshingCell&, ResultingCell&) const { }
		} cellConverter;
		
		// try to repeat Triangulation_2 copy constructor as much as possible
		result.set_infinite_vertex(result.tds().copy_tds(
				mt.tds(), mt.infinite_vertex(), vertexConverter, cellConverter));
	}
	
	
};


}

#endif // LIBCGALMESH_CGAL2DMESHER_HPP
