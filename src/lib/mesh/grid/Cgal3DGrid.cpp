#include <lib/mesh/grid/Cgal3DGrid.hpp>
#include <lib/util/FileUtils.hpp>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>


using namespace gcm;


Cgal3DGrid::
Cgal3DGrid(const Task& task) :
	UnstructuredGrid(task),
	effectiveSpatialStep(task.cgal3DGrid.spatialStep), 
	movable(task.cgal3DGrid.movable) {

	triangulate(task.cgal3DGrid);
	
	vertexHandles.resize(triangulation.number_of_vertices());
	size_t vertexIndex = 0;
	for (auto it  = triangulation.finite_vertices_begin();
	          it != triangulation.finite_vertices_end(); it++) {
		it->info() = vertexIndex;
		vertexHandles[vertexIndex] = it;
		vertexIndex++;
	}

	markInnersAndBorders();
}


std::vector<Cgal3DGrid::Iterator> Cgal3DGrid::
findNeighborVertices(const Iterator& it) const {	

	std::vector<VertexHandle> neighborVertices;
	triangulation.finite_adjacent_vertices(vertexHandles[getIndex(it)],
			std::back_inserter(neighborVertices));

	std::vector<Iterator> ans;
	for (const auto vertex : neighborVertices) {
		ans.push_back(getIterator(vertex));
	}
	return ans;
}


Cgal3DGrid::Cell Cgal3DGrid::
findOwnerCell(const Iterator& /*it*/, const Real3& /*shift*/) const {
	THROW_UNSUPPORTED("TODO");
}


Cgal3DGrid::Cell Cgal3DGrid::
locateOwnerCell(const Iterator& it, const Real3& shift) const {
	VertexHandle beginVertex = vertexHandles[getIndex(it)];
	CgalPoint3 query = beginVertex->point() + cgalVector3(shift);
	CellHandle ownerCell = triangulation.locate(query, beginVertex->cell());
	return createTetrahedron(ownerCell);
}


Cgal3DGrid::Face Cgal3DGrid::
findCrossingBorder(const Iterator& start, const Real3& direction) const {
	VertexHandle v = vertexHandles[getIndex(start)];
	CellHandle startCell = findCrossedCell(v, direction);
	return findCrossingBorder(v, startCell, direction);
}


std::vector<Cgal3DGrid::Iterator> Cgal3DGrid::
findBorderNeighbors(const Iterator& it) const {
	/// only for border nodes
	assert_true(isBorder(it));	
	THROW_UNSUPPORTED("TODO");
}


Cgal3DGrid::Iterator Cgal3DGrid::
findVertexByCoordinates(const Real3& coordinates) const {
	for (const auto& it : *this) {
		if (coords(it) == coordinates) {
			return it;
		}
	}
	THROW_INVALID_ARG("There isn't a vertex with such coordinates");
}


void Cgal3DGrid::
triangulate(const Task::Cgal3DGrid& task) {

	// Types connected with meshing
	// Initial mesh "domain" - from what CGAL builds the triangulation
	typedef CGAL::Polyhedron_3<K>                                         Polyhedron;
	typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K>                 PolyhedralMeshDomain;
	// Triangulation used by mesher
	typedef CGAL::Mesh_triangulation_3<PolyhedralMeshDomain>::type        MeshingTriangulation;
	typedef CGAL::Mesh_complex_3_in_triangulation_3<MeshingTriangulation> C3t3;
	typedef typename C3t3::Subdomain_index                                SubdomainIndex;
	typedef typename MeshingTriangulation::Vertex                         MeshingVertex;
	typedef typename MeshingTriangulation::Cell                           MeshingCell;
	// Meshing criteria
	typedef CGAL::Mesh_criteria_3<MeshingTriangulation>                   MeshingCriteria;

	Polyhedron polyhedron;
	FileUtils::readFromTextFile(task.polyhedronFileName, polyhedron);
	assert_true(polyhedron.is_valid());
	
	// Create initial mesh domain
	PolyhedralMeshDomain domain(polyhedron);
	
	MeshingCriteria meshingCriteria(
			CGAL::parameters::facet_angle=25,
			CGAL::parameters::cell_radius_edge_ratio=3,
			CGAL::parameters::facet_size=task.spatialStep,
			CGAL::parameters::facet_distance=task.spatialStep / 4,
			CGAL::parameters::cell_size=task.spatialStep);

	LOG_DEBUG("Meshing the triangulation...");
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, meshingCriteria,
			CGAL::parameters::no_perturb(), CGAL::parameters::no_exude());
	MeshingTriangulation mt = c3t3.triangulation();
	assert_true(mt.is_valid());
	
	// load triangulation from mesher
	// converter structures for CGAL copy_tds function
	struct VertexConverter {
		VertexT operator()(const MeshingVertex& src) const {
			return VertexT(src.point());
		}
		void operator()(const MeshingVertex&, VertexT&) const { }
	} vertexConverter;
	struct CellConverter {
		CellT operator()(const MeshingCell& orig) const {
			CellT cellT;
			bool cellIsInComplex = !(orig.subdomain_index() == SubdomainIndex());
			cellT.info() = cellIsInComplex ? 1 : 0;
			return cellT;
		}
		void operator()(const MeshingCell&, CellT&) const { }
	} cellConverter;
	
	// try to repeat Triangulation_3 copy constructor as much as possible
	triangulation.set_lock_data_structure(mt.get_lock_data_structure());
	triangulation.set_infinite_vertex(triangulation.tds().copy_tds(mt.tds(),
			mt.infinite_vertex(), vertexConverter, cellConverter));
	
	assert_true(triangulation.is_valid());
	LOG_DEBUG("Number of vertices after meshing: " << triangulation.number_of_vertices());
	LOG_DEBUG("Number of cells after meshing: " << triangulation.number_of_cells());
}


void Cgal3DGrid::markInnersAndBorders() {
	/// insert indices of border vertices into borderIndices
	/// and indices of inner vertices into innerIndices
	borderIndices.clear();
	innerIndices.clear();
	
	for (auto v  = triangulation.finite_vertices_begin();
	          v != triangulation.finite_vertices_end(); ++v) {

		bool isBorderNode = false;
		std::vector<CellHandle> incidentCells;
		triangulation.incident_cells(v, std::back_inserter(incidentCells));
		for (const auto cell : incidentCells) {
			if ( !isInDomain(cell) ) {
				isBorderNode = true;
				break;
			}
		}
		
		if (isBorderNode) {
			borderIndices.insert(getIndex(getIterator(v)));
		} else {
			innerIndices.insert(getIndex(getIterator(v)));
		}
		
	}
	
	assert_eq(borderIndices.size() + innerIndices.size(), sizeOfAllNodes());
	LOG_DEBUG("Number of border vertices: " << borderIndices.size());
	LOG_DEBUG("Number of inner vertices: " << innerIndices.size());
}


Cgal3DGrid::CellHandle Cgal3DGrid::
findCrossedCell(const VertexHandle start, const Real3& direction) const {
	/// Choose among incident cells of the given point that one which is
	/// crossed by the line from that point in specified direction.

	std::vector<CellHandle> finiteIncidentCells;
	triangulation.finite_incident_cells(start, std::back_inserter(finiteIncidentCells));
	for (int n = 1; n < 20; n++) {
		CgalPoint3 q = start->point() + cgalVector3(direction / pow(2, n));
		for (const auto cell : finiteIncidentCells) {
			if (!triangulation.tetrahedron(cell).has_on_unbounded_side(q)) {
				return cell;
			}
		}
	}
	THROW_BAD_MESH("Exceed number of iterations in findCrossedCell");
}


std::vector<Cgal3DGrid::VertexHandle> Cgal3DGrid::
commonVertices(const CellHandle& a, const CellHandle& b) const {
	/// return common vertices of given cells
	std::vector<VertexHandle> ans;
	for (int i = 0; i < 4; i++) {
		VertexHandle candidate = a->vertex(i);
		if (b->has_vertex(candidate)) {
			ans.push_back(candidate);
		}
	}
	return ans;
}


Cgal3DGrid::Face Cgal3DGrid::
findCrossingBorder(const VertexHandle& /*v*/, const CellHandle& /*f*/,
		const Real3& /*direction*/) const {
	/// given a start point and direction, go along this line until
	/// intersection with border;
	/// @note cell is incident for given vertex,
	/// given line must go across given cell
	/// @return found border as pair of its vertices
	
	THROW_UNSUPPORTED("TODO");
}


void Cgal3DGrid::
printCell(const CellHandle& f, const std::string& name) const {
	/// debugging helper
	SUPPRESS_WUNUSED(name);
	LOG_DEBUG("Cell " << name << ":");
	for (int i = 0; i < 4; i++) {
		if (triangulation.is_infinite(f->vertex(i))) {
			LOG_DEBUG("\nINFINITE\n");
		} else {
			LOG_DEBUG(real3(f->vertex(i)->point()));
		}
	}
}








