#ifndef LIBGCM_GRIDCHARACTERISTICMETHODCUBICGRID_HPP
#define LIBGCM_GRIDCHARACTERISTICMETHODCUBICGRID_HPP

#include <lib/numeric/gcm/GridCharacteristicMethod.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>


namespace gcm {

/**
 * Grid-characteristic method specialization for meshes based on CubicGrid
 */
template<typename TModel, typename TMaterial, int Dimensionality>
class GridCharacteristicMethod<TModel, CubicGrid<Dimensionality>, TMaterial> {
public:
	typedef CubicGrid<Dimensionality>            Grid;
	typedef DefaultMesh<TModel, Grid, TMaterial> Mesh;
	typedef typename Mesh::Matrix                Matrix;
	typedef typename Mesh::PdeVector             PdeVector;
	typedef typename Mesh::Iterator              Iterator;

	/**
	 * Do grid-characteristic stage of splitting method
	 * @param s direction aka stage
	 * @param timeStep time step
	 * @param mesh mesh to perform calculation
	 */
	void stage(const int s, const real& timeStep, Mesh& mesh) const {
		for (auto it : mesh) {
			mesh._pdeNew(it) = localGcmStep(
					mesh.matrices(it)->m[s].U1,
					mesh.matrices(it)->m[s].U,
					interpolateValuesAround(mesh, s, it,
							crossingPoints(it, s, timeStep, mesh)));
		}
	}

public: // TODO - private (rewrite test)
	/** Points where characteristics from next time layer cross current time layer */
	PdeVector crossingPoints(const Iterator& it, const int s,
	                         const real& timeStep, const Mesh& mesh) const {
		return -timeStep * linal::diag(mesh.matrices(it)->m[s].L);
	}
	
	
	/**
	 * Interpolate nodal values in specified points.
	 * Interpolated value for k-th point in vector dx
	 * stored in k-th column of returned Matrix.
	 * @param mesh mesh to perform interpolation on
	 * @param s direction
	 * @param it index-iterator of node
	 * @param dx Vector of distances from reference node on which
	 * values should be interpolated
	 * @return Matrix with interpolated nodal values in columns
	 */
	Matrix interpolateValuesAround(const Mesh& mesh, const int s,
	                               const Iterator& it, const PdeVector& dx) const {
		Matrix ans;
		std::vector<PdeVector> src( (size_t) (mesh.borderSize + 1) );
		Iterator shift = Iterator::Zeros();
		for (int k = 0; k < PdeVector::M; k++) {
			shift(s) = (dx(k) > 0) ? 1 : -1;
			for (int i = 0; i < src.size(); i++) {
				src[(size_t)i] = mesh.pde(it + shift * i);
			}
			ans.setColumn(k, EqualDistanceLineInterpolator<PdeVector>::
					minMaxInterpolate(src, fabs(dx(k)) / mesh.h(s)));
		}
		return ans;
	}
	
};


}

#endif // LIBGCM_GRIDCHARACTERISTICMETHODCUBICGRID_HPP
