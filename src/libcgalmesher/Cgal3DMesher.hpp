#ifndef LIBCGALMESH_CGAL3DMESHER_HPP
#define LIBCGALMESH_CGAL3DMESHER_HPP

#include <string>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>



namespace cgalmesher {


/**
 * 3D mesher by CGAL library.
 * The reason why each meshing function is in separate cpp file --
 * horrible memory usage (gcc-4,-5) and low speed of compilation 
 * because of swapping on even 4Gb system.
 * (Headers only, boost, include all to each other ..)
 * Including headers of CGAL meshers to cpp only,
 * we can compile a separate shared library
 * just once and forget about it developing main project.
 */
class Cgal3DMesher {
public:
	
	typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
	typedef CGAL::Triangulation_vertex_base_with_info_3<size_t, K> Vb;
	typedef CGAL::Triangulation_cell_base_with_info_3<size_t, K>   Cb;
	typedef CGAL::Triangulation_data_structure_3<Vb, Cb>           Tds;
	
	/// The result of meshing to be copied to some resulting triangulation
	typedef CGAL::Delaunay_triangulation_3<K, Tds> IntermediateTriangulation;
	typedef IntermediateTriangulation::Vertex      Vertex;
	typedef IntermediateTriangulation::Cell        Cell;
	
	
	/// @name Converter structures for CGAL copy_tds function --
	/// convertion between vertices and cells from
	/// different types of triangulations
	/// @{
	
	template<typename Src, typename Res>
	class DefaultVertexConverter {
	public:
		Res operator()(const Src& src) const {
			return Res(src.point()); // copy coordinates
		}
		void operator()(const Src&, Res&) const { }
	};
	
	template<typename Src, typename Res>
	struct DefaultCellConverter {
		Res operator()(const Src& src) const {
			Res res;
			// copy information about empty/non-empty space in the cell
			res.info().setGridId(src.info());
			return res;
		}
		void operator()(const Src&, Res&) const { }
	};
	
	template<typename Src, typename Res>
	struct SubdomainCellConverter {
		typedef typename Src::Subdomain_index SubdomainIndex;
		
		Res operator()(const Src& src) const {
			Res res;
			res.info() = (size_t)(-1);
			// copy information about empty/non-empty space in the cell
			if (src.subdomain_index() != SubdomainIndex()) {
				res.info() = 0;
			}
			return res;
		}
		void operator()(const Src&, Res&) const { }
	};
	
	/// @}
	
	
	/**
	 * Build the grid on given geometry
	 * @param spatialStep         effective spatial step
	 * @param detectSharpEdges    use true for figures with sharp edges
	 * @param polyhedronFileName  file with polyhedron to construct the grid from
	 * @param result              triangulation to write the result in
	 * @tparam ResultingTriangulation type of the triangulation to write result in
	 * @tparam CellConverter          see DefaultCellConverter
	 * @tparam VertexConverter        see DefaultVertexConverter
	 */
	template<
			typename ResultingTriangulation,
			template<typename, typename> class CellConverter = DefaultCellConverter,
			template<typename, typename> class VertexConverter = DefaultVertexConverter
			>
	static void triangulate(const double spatialStep, const bool detectSharpEdges,
			const std::string polyhedronFileName, ResultingTriangulation& result) {
		
		IntermediateTriangulation intermediateTriangulation;
		
		if (detectSharpEdges) {
			intermediateTriangulation = triangulateWithEdges(polyhedronFileName, spatialStep);
		} else {
			intermediateTriangulation = triangulateWithoutEdges(polyhedronFileName, spatialStep);
		}
		
		copyTriangulation<IntermediateTriangulation, ResultingTriangulation,
				CellConverter, VertexConverter>(intermediateTriangulation, result);
	}
	
	
private:
	
	/**
	 * Triangulation process that meshes cube figure into the cube figure and 
	 * tetrahedron into tetrahedron one - initial edges of the figure stay preserved.
	 * However, resulting triangulation may not be Delaunay triangulation -
	 * property of empty sphere can be violated on sharp edges.
	 * @param polyhedronFileName  file with polyhedron to construct the grid from
	 * @param spatialStep         effective spatial step
	 */
	static IntermediateTriangulation triangulateWithEdges(
			const std::string polyhedronFileName, const double spatialStep);
	
	
	/**
	 * Triangulation process that smoothes any sharp edges of initial figure.
	 * However, resulting triangulation is strictly Delaunay triangulation.
	 * @param polyhedronFileName  file with polyhedron to construct the grid from
	 * @param spatialStep         effective spatial step
	 */
	static IntermediateTriangulation triangulateWithoutEdges(
			const std::string polyhedronFileName, const double spatialStep);
	
	
	/**
	 * Copy CGAL triangulations of different types
	 */
	template<
			typename InputTriangulation,
			typename OutputTriangulation,
			template<typename, typename> class CellConverter,
			template<typename, typename> class VertexConverter
			>
	static void copyTriangulation(const InputTriangulation& input,
			OutputTriangulation& output) {
		
		typedef typename InputTriangulation::Cell     InputCell;
		typedef typename InputTriangulation::Vertex   InputVertex;
		typedef typename OutputTriangulation::Cell    OutputCell;
		typedef typename OutputTriangulation::Vertex  OutputVertex;
		
		CellConverter<InputCell, OutputCell>        cellConverter;
		VertexConverter<InputVertex, OutputVertex>  vertexConverter;
		
		copyTriangulation(input, output, cellConverter, vertexConverter);
	}
	
	
	/**
	 * Copy CGAL triangulations of different types
	 */
	template<
			typename InputTriangulation,
			typename OutputTriangulation,
			typename CellConverter,
			typename VertexConverter
			>
	static void copyTriangulation(
			const InputTriangulation& input, OutputTriangulation& output,
			const CellConverter& cellConverter,
			const VertexConverter& vertexConverter) {
		// try to repeat Triangulation_3 copy constructor as much as possible
		output.set_lock_data_structure(input.get_lock_data_structure());
		output.set_infinite_vertex(output.tds().copy_tds(
				input.tds(), input.infinite_vertex(), vertexConverter, cellConverter));
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
