#ifndef LIBGCM_GCMHANDLER_HPP
#define LIBGCM_GCMHANDLER_HPP

#include <lib/linal/linal.hpp>
#include <lib/mesh/DefaultMesh.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>
#include <lib/mesh/grid/Cgal2DGrid.hpp>
#include <lib/numeric/interpolation/interpolation.hpp>


namespace gcm {
	/**
	 * Bridge between grids of different types and grid-characteristic method
	 */
	template<typename TModel, typename TGrid>
	struct GcmHandler {
		typedef DefaultMesh<TModel, TGrid>           Mesh;
		typedef typename Mesh::Matrix                Matrix;
		typedef typename Mesh::PdeVector             PdeVector;
		typedef typename Mesh::Iterator              Iterator;

		/**
		 * Interpolate nodal values in specified points.
		 * Interpolated value for k-th point in vector dx are
		 * stored in k-th column of returned Matrix.
		 * @param stage direction
		 * @param it index-iterator of node
		 * @param dx Vector of distances from reference node on which
		 * values should be interpolated
		 * @return Matrix with interpolated nodal values in columns
		 */
		static Matrix interpolateValuesAround(const Mesh& mesh, const int stage, const Iterator& it,
		                                      const PdeVector& dx);
	};

	template<typename TModel>
	struct GcmHandler<TModel, CubicGrid> {
		typedef DefaultMesh<TModel, CubicGrid>       Mesh;
		typedef typename Mesh::Matrix                Matrix;
		typedef typename Mesh::PdeVector             PdeVector;
		typedef typename Mesh::Iterator              Iterator;

		/** See the comment under */
		static Matrix interpolateValuesAround(const Mesh& mesh, const int stage,
		                                      const Iterator& it, const PdeVector& dx) {
			Matrix ans;
			std::vector<PdeVector> src( (size_t) (mesh.getAccuracyOrder() + 1) );
			PdeVector res;
			Iterator shift({0, 0, 0}, mesh.getSizes());
			for (int k = 0; k < PdeVector::M; k++) {
				shift(stage) = (dx(k) > 0) ? 1 : -1;
				for (int i = 0; i < src.size(); i++) {
					src[(size_t)i] = mesh.pde(it + shift * i);
				}
				EqualDistanceLineInterpolator<PdeVector>::minMaxInterpolate(res, src, fabs(dx(k)) / mesh.h(stage));
				ans.setColumn(k, res);
			}
			return ans;
		};
	};

	template<typename TModel>
	struct GcmHandler<TModel, Cgal2DGrid> {
		typedef DefaultMesh<TModel, Cgal2DGrid>      Mesh;
		typedef typename Mesh::Matrix                Matrix;
		typedef typename Mesh::PdeVector             PdeVector;
		typedef typename Mesh::Iterator              Iterator;
		typedef typename Mesh::Triangle              Triangle;

		/** See the comment under */
		static Matrix interpolateValuesAround(const Mesh& mesh, const int stage,
		                                      const Iterator& it, const PdeVector& dx) {
			Matrix ans;
			for (int k = 0; k < PdeVector::M; k++) {
				linal::Vector2 shift({(stage == 0) * dx(k), (stage == 1) * dx(k)});
				auto t = mesh.findOwnerTriangle(it, shift);
				PdeVector u; linal::clear(u);
				if (t.inner) {
					u = TriangleInterpolator<PdeVector>::interpolate(
							mesh.coords2d(t.p[0]), mesh.pde(t.p[0]),
							mesh.coords2d(t.p[1]), mesh.pde(t.p[1]),
							mesh.coords2d(t.p[2]), mesh.pde(t.p[2]),
							mesh.coords2d(it) + shift);
				}
				ans.setColumn(k, u);
			}
			return ans;
		};
	};
}

#endif // LIBGCM_GCMHANDLER_HPP
