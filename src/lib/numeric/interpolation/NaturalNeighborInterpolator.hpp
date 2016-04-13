#ifndef LIBGCM_NATURALNEIGHBORINTERPOLATOR_HPP
#define LIBGCM_NATURALNEIGHBORINTERPOLATOR_HPP

#include <CGAL/Kernel/global_functions.h>

#include <list>

#include <lib/linal/linal.hpp>

namespace gcm {

/**
 * Implements "Natural neighbor" interpolation method.
 * Most code are copypasted from CGAL/natural_neighbor_coordinates_2.h
 * and modified for our purposes.
 * This type of interpolation is much slower than usual, but smoother.
 */
template<typename TMesh>
class NaturalNeighborInterpolator {
public:
	typedef typename TMesh::Grid                Grid;
	typedef typename Grid::CDT                  CDT;
	typedef typename Grid::K                    Kernel;
	typedef typename Kernel::FT                 Coordinate;
	typedef typename CDT::Vertex_handle         VertexHandle;
	typedef typename CDT::Face_handle           FaceHandle;
	typedef typename CDT::Point                 Point;
	typedef typename CDT::Edge                  Edge;

	typedef std::pair<VertexHandle, Coordinate> NaturalNeighbor;
	typedef std::vector<NaturalNeighbor>        NaturalNeighbors;

	typedef typename TMesh::PdeVector           Value;


	static Value firstOrderSibsonInterpolate(const TMesh& mesh, const Real2& query,
	                                         const VertexHandle& nearVertex) {
		NaturalNeighbors neighbors = naturalNeighborCoordinates(
		        mesh.getTriangulation(), Point(query(0), query(1)), nearVertex);

		/// first order classic Sibson interpolation
		Value ans = Value::zeros();
		for (const auto& neighbor : neighbors) {
			const Value neighborValue =
			        mesh.pde(mesh.getIterator(neighbor.first));
			ans += neighborValue * neighbor.second;
		}

		return ans;
	}

	/** @return normalized natural neighbor coordinates for query */
	static NaturalNeighbors
	naturalNeighborCoordinates(const CDT& triangulation, const Point& query,
	                           const VertexHandle& nearVertex) {

		NaturalNeighbors neighbors;
		Coordinate norm = naturalNeighborCoordinates(triangulation, query,
				neighbors, nearVertex->incident_faces());

		for (auto& n : neighbors) {
			n.second /= norm;
		}
		return neighbors;
	}


private:
	static Coordinate
	naturalNeighborCoordinates(const CDT& triangulation, const Point p, 
	                           NaturalNeighbors& neighbors, const FaceHandle start) {

		typename CDT::Locate_type lt;
		int li;
		FaceHandle fh = triangulation.locate(p, lt, li, start);

		assert_false(lt == CDT::OUTSIDE_AFFINE_HULL || lt == CDT::OUTSIDE_CONVEX_HULL);

		if ((lt == CDT::EDGE &&
		     (triangulation.is_infinite(fh) ||
		      triangulation.is_infinite(fh->neighbor(li))))) {
			VertexHandle v1 = fh->vertex(triangulation.cw(li));
			VertexHandle v2 = fh->vertex(triangulation.ccw(li));

			Point p1(v1->point()), p2(v2->point());

			Coordinate coef1(0);
			Coordinate coef2(0);
			typename CDT::Geom_traits::Equal_x_2 equal_x_2;
			if (!equal_x_2(p1, p2)) {
				coef1 = (p.x() - p2.x()) / (p1.x() - p2.x());
				coef2 = 1 - coef1;
				neighbors = {{v1, coef1}, {v2, coef2}};
			} else {
				coef1 = (p.y() - p2.y()) / (p1.y() - p2.y());
				coef2 = 1 - coef1;
				neighbors = {{v1, coef1}, {v2, coef2}};
			}

			return coef1 + coef2;
		}

		if (lt == CDT::VERTEX) {
			neighbors = {{fh->vertex(li), Coordinate(1)}};
			return Coordinate(1);
		}

		std::list<Edge> hole;
		triangulation.get_boundary_of_conflicts(p, std::back_inserter(hole), fh);
		return naturalNeighborCoordinates(triangulation, p, neighbors, hole);
	}


	static Coordinate
	naturalNeighborCoordinates(const CDT& triangulation, const Point& p,
	                           NaturalNeighbors& neighbors, const std::list<Edge>& hole) {
		/// function call if the conflict zone is known
		std::vector<Point> vor(3);
		Coordinate area_sum(0);

		auto hit = hole.end();
		--hit;
		// in the beginning: prev is the "last" vertex of the hole:
		// later: prev is the last vertex processed (previously)
		VertexHandle prev = hit->first->vertex(triangulation.cw(hit->second));
		hit = hole.begin();

		while (hit != hole.end()) {
			Coordinate area(0);
			VertexHandle current = hit->first->vertex(triangulation.cw(hit->second));

			vor[0] = triangulation.geom_traits().construct_circumcenter_2_object()
					(current->point(),
					 hit->first->vertex(triangulation.ccw(hit->second))->point(), p);

			auto fc = triangulation.incident_faces(current, hit->first);
			++fc;
			vor[1] = CGAL::circumcenter(triangulation.triangle(fc));

			while (!fc->has_vertex(prev)) {
				++fc;
				vor[2] = CGAL::circumcenter(triangulation.triangle(fc));

				area += polygon_area_2(vor.begin(),
				                       vor.end(), triangulation.geom_traits());

				vor[1] = vor[2];
			}
			vor[2] = triangulation.geom_traits().construct_circumcenter_2_object() 
					(prev->point(), current->point(), p);
			area += polygon_area_2(vor.begin(), vor.end(), triangulation.geom_traits());


			neighbors.push_back({current, area});
			area_sum += area;

			// update prev and hit:
			prev = current;
			++hit;
		}
		return area_sum;
	}
	
	
};


}

#endif // LIBGCM_NATURALNEIGHBORINTERPOLATOR_HPP
