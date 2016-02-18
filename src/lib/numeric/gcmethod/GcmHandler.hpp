#ifndef LIBGCM_GCMHANDLER_HPP
#define LIBGCM_GCMHANDLER_HPP

#include <lib/linal/linal.hpp>
#include <lib/mesh/DefaultMesh.hpp>
#include <lib/mesh/CubicGrid.hpp>
#include <lib/mesh/Cgal2DGrid.hpp>
#include <lib/numeric/interpolation/EqualDistanceLineInterpolator.hpp>


namespace gcm {
	/**
	 * Bridge between grids of different types and grid-characteristic method
	 */
	template<typename TModel, typename TGrid>
	struct GcmHandler {
		typedef DefaultMesh<TModel, TGrid>           Grid;
		typedef typename Grid::Matrix                Matrix;
		typedef typename Grid::PdeVector             PdeVector;
		typedef typename Grid::Iterator              Iterator;

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
		static Matrix interpolateValuesAround(const Grid& grid, const int stage, const Iterator& it,
		                                      const PdeVector& dx);
	};

	template<typename TModel>
	struct GcmHandler<TModel, CubicGrid> {
		typedef DefaultMesh<TModel, CubicGrid>       Grid;
		typedef typename Grid::Matrix                Matrix;
		typedef typename Grid::PdeVector             PdeVector;
		typedef typename Grid::Iterator              Iterator;

		/** See the comment under */
		static Matrix interpolateValuesAround(const Grid& grid, const int stage,
		                                      const Iterator& it, const PdeVector& dx) {
			Matrix ans;
			std::vector<PdeVector> src( (size_t) (grid.getAccuracyOrder() + 1) );
			PdeVector res;
			Iterator shift({0, 0, 0}, grid.getSizes());
			for (int k = 0; k < PdeVector::M; k++) {
				shift(stage) = (dx(k) > 0) ? 1 : -1;
				for (int i = 0; i < src.size(); i++) {
					src[(size_t)i] = grid.pde(it + shift * i);
				}
				EqualDistanceLineInterpolator<PdeVector>::minMaxInterpolate(res, src, fabs(dx(k)) / grid.h(stage));
				ans.setColumn(k, res);
			}
			return ans;
		};
	};

	template<typename TModel>
	struct GcmHandler<TModel, Cgal2DGrid> {
		typedef DefaultMesh<TModel, Cgal2DGrid>      Grid;
		typedef typename Grid::Matrix                Matrix;
		typedef typename Grid::PdeVector             PdeVector;
		typedef typename Grid::Iterator              Iterator;

		/** See the comment under */
		static Matrix interpolateValuesAround(const Grid& grid, const int stage,
		                                      const Iterator& it, const PdeVector& dx) {
			SUPPRESS_WUNUSED(grid);
			SUPPRESS_WUNUSED(stage);
			SUPPRESS_WUNUSED(it);
			SUPPRESS_WUNUSED(dx);
			Matrix ans;
			linal::clear(ans);
			return ans;
		};
	};
}

#endif // LIBGCM_GCMHANDLER_HPP
