#include <libcgalmesh/Cgal3DMesher.hpp>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

using namespace cgalmesh;


void Cgal3DMesher::
triangulateWithoutEdges(const std::string polyhedronFileName,
			const double spatialStep, ResultingTriangulation& result) {

	// Types connected with meshing
	// Initial mesh "domain" - from what CGAL builds the triangulation
	typedef CGAL::Polyhedron_3<K>                                         Polyhedron;
	typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K>                 MeshDomain;
	// Triangulation used by mesher
	typedef CGAL::Mesh_triangulation_3<MeshDomain>::type                  MeshingTriangulation;
	typedef CGAL::Mesh_complex_3_in_triangulation_3<MeshingTriangulation> C3t3;
	// Meshing criteria
	typedef CGAL::Mesh_criteria_3<MeshingTriangulation>                   MeshingCriteria;

	Polyhedron polyhedron;
	readFromTextFile(polyhedronFileName, polyhedron);
	assert(polyhedron.is_valid());
			
	// Create initial mesh domain
	MeshDomain domain(polyhedron);
	
	MeshingCriteria meshingCriteria(
			CGAL::parameters::facet_topology = CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH,
			CGAL::parameters::facet_angle = 25,
			CGAL::parameters::cell_radius_edge_ratio = 3,
			CGAL::parameters::facet_size = spatialStep,
			CGAL::parameters::facet_distance = spatialStep / 4,
			CGAL::parameters::cell_size = spatialStep);

	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, meshingCriteria,
			CGAL::parameters::no_perturb(), CGAL::parameters::no_exude());
	MeshingTriangulation mt = c3t3.triangulation();
	assert(mt.is_valid());
	
	copyTriangulation(mt, result);
	assert(result.is_valid());
}

