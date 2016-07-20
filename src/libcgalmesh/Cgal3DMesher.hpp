#ifndef LIBCGALMESH_CGAL3DMESHER_HPP
#define LIBCGALMESH_CGAL3DMESHER_HPP

#include <string>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>


namespace cgalmesh {

/**
 * 3D mesher by CGAL library
 */
class Cgal3DMesher {
public:
	typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
	typedef CGAL::Triangulation_vertex_base_with_info_3<size_t, K> Vb;
	typedef CGAL::Triangulation_cell_base_with_info_3<size_t, K>   Cb;
	typedef CGAL::Triangulation_data_structure_3<Vb, Cb>           Tds;
	typedef CGAL::Delaunay_triangulation_3<K, Tds>                 ResultingTriangulation;
		
	typedef ResultingTriangulation::Vertex                         ResultingVertex;
	typedef ResultingTriangulation::Cell                           ResultingCell;

	/**
	 * Build the grid on given geometry
	 * @param spatialStep         effective spatial step
	 * @param detectSharpEdges    use true for figures with sharp edges
	 * @param polyhedronFileName  file with polyhedron to construct the grid from
	 * @param result              triangulation to write the result in
	 */
	static void triangulate(const double spatialStep, const bool detectSharpEdges,
			const std::string polyhedronFileName, ResultingTriangulation& result) {
		
		if (detectSharpEdges) {
			triangulateWithEdges(polyhedronFileName, spatialStep, result);
		} else {
			triangulateWithoutEdges(polyhedronFileName, spatialStep, result);
		}
	}

	/**
	 * Triangulation process that meshes cube figure into the cube figure and 
	 * tetrahedron into tetrahedron one - initial edges of the figure stay preserved.
	 * However, resulting triangulation may not be Delaunay triangulation -
	 * property of empty sphere can be violated on sharp edges.
	 * @param polyhedronFileName  file with polyhedron to construct the grid from
	 * @param spatialStep         effective spatial step
	 * @param result              triangulation to write the result in
	 */
	static void triangulateWithEdges(const std::string polyhedronFileName,
			const double spatialStep, ResultingTriangulation& result);

	/**
	 * Triangulation process that smoothes any sharp edges of initial figure.
	 * However, resulting triangulation is strictly Delaunay triangulation.
	 * @param polyhedronFileName  file with polyhedron to construct the grid from
	 * @param spatialStep         effective spatial step
	 * @param result              triangulation to write the result in
	 */
	static void triangulateWithoutEdges(const std::string polyhedronFileName,
			const double spatialStep, ResultingTriangulation& result);
	

private:

	/** Load result triangulation from mesher triangulation */
	template<typename MeshingTriangulation>
	static void copyTriangulation(const MeshingTriangulation& mt, 
			ResultingTriangulation& result) {
		
		typedef typename MeshingTriangulation::Vertex    MeshingVertex;
		typedef typename MeshingTriangulation::Cell      MeshingCell;
		typedef typename MeshingCell::Subdomain_index    SubdomainIndex;
		
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
				bool cellIsInComplex = !(src.subdomain_index() == SubdomainIndex());
				cellT.info() = cellIsInComplex ? 1 : 0;
				return cellT;
			}
			void operator()(const MeshingCell&, ResultingCell&) const { }
		} cellConverter;
		
		// try to repeat Triangulation_3 copy constructor as much as possible
		result.set_lock_data_structure(mt.get_lock_data_structure());
		result.set_infinite_vertex(result.tds().copy_tds(
				mt.tds(), mt.infinite_vertex(), vertexConverter, cellConverter));
	}
	
	
	template<typename PlaceToWriteInput>
	static void readFromTextFile(const std::string& fileName,
			PlaceToWriteInput& placeToWriteInput) {
		
		std::ifstream inputFileStream(fileName);
		assert(inputFileStream.is_open());
		
		inputFileStream >> placeToWriteInput;
		assert(inputFileStream.good());
		inputFileStream.close();
	}
	
};


}

#endif // LIBCGALMESH_CGAL3DMESHER_HPP
