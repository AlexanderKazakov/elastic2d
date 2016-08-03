#include <libcgalmesher/Cgal3DMesher.hpp>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>


using namespace cgalmesher;


Cgal3DMesher::IntermediateTriangulation
Cgal3DMesher::
triangulateWithEdges(const std::string polyhedronFileName,
			const double spatialStep) {
	
	typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> MeshDomainWithEdges;
	typedef CGAL::Mesh_triangulation_3<
			MeshDomainWithEdges>::type MeshTriangulationWithEdges;
	typedef CGAL::Mesh_complex_3_in_triangulation_3<
			MeshTriangulationWithEdges,
			MeshDomainWithEdges::Corner_index,
			MeshDomainWithEdges::Curve_segment_index> C3t3;
	typedef CGAL::Mesh_criteria_3<MeshTriangulationWithEdges> MeshingCriteria;
	
	// Create initial mesh domain
	MeshDomainWithEdges domain(polyhedronFileName.c_str());
	// Identify sharp edges (features)
	domain.detect_features();
	
	MeshingCriteria meshingCriteria(
			CGAL::parameters::edge_size = spatialStep,
			CGAL::parameters::facet_topology = CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH,
			CGAL::parameters::facet_angle = 25,
			CGAL::parameters::cell_radius_edge_ratio = 3,
			CGAL::parameters::facet_size = spatialStep,
			CGAL::parameters::facet_distance = spatialStep / 4,
			CGAL::parameters::cell_size = spatialStep);
	
	// meshing
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, meshingCriteria,
			CGAL::parameters::no_perturb(), CGAL::parameters::no_exude());
	
	MeshTriangulationWithEdges mt = c3t3.triangulation();
	assert(mt.is_valid());
	
	IntermediateTriangulation result;
	copyTriangulation<
			MeshTriangulationWithEdges, IntermediateTriangulation,
			SubdomainCellConverter, DefaultVertexConverter>(
					mt, result);
	
	return result;
}

