#ifndef LIBGCM_CUBIC_GRIDCHARACTERISTICMETHOD_HPP
#define LIBGCM_CUBIC_GRIDCHARACTERISTICMETHOD_HPP

#include <libgcm/grid/AbstractGrid.hpp>
#include <libgcm/util/math/GridCharacteristicMethod.hpp>
#include <libgcm/util/math/interpolation/interpolation.hpp>


namespace gcm {
namespace cubic {


class GridCharacteristicMethodBase {
public:
	virtual void stage(
			const int s, const real& timeStep, AbstractGrid& mesh_) const = 0;
};



/**
 * Grid-characteristic method for meshes based on CubicGrid
 */
template<typename Mesh>
class GridCharacteristicMethod : public GridCharacteristicMethodBase {
public:
	typedef typename Mesh::Matrix                Matrix;
	typedef typename Mesh::PdeVector             PdeVector;
	typedef typename Mesh::Iterator              Iterator;
	
	
	GridCharacteristicMethod(const Task&) { }
	
	
	/**
	 * Do grid-characteristic stage of splitting method
	 * @param s stage (here is equal to direction 0 == X, 1 == Y, 2 == Z)
	 * @param timeStep time step
	 * @param mesh mesh to perform calculation
	 * collinear with coordinate axes are available
	 */
	virtual void stage(
			const int s, const real& timeStep, AbstractGrid& mesh_) const override {
		Mesh& mesh = dynamic_cast<Mesh&>(mesh_);
		for (auto it : mesh) {
			mesh._pdeNew(s, it) = localGcmStep(
					mesh.matrices(it)->m[s].U1,
					mesh.matrices(it)->m[s].U,
					interpolateValuesAround(mesh, s, it,
							crossingPoints(it, s, timeStep, mesh)));
		}
	}
	
	
	/** Points where characteristics from next time layer cross current time layer */
	static PdeVector crossingPoints(const Iterator& it, const int s,
			const real& timeStep, const Mesh& mesh) {
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
	static Matrix interpolateValuesAround(const Mesh& mesh, const int s,
			const Iterator& it, const PdeVector& dx) {
		Matrix ans;
		std::vector<PdeVector> src( (size_t) (mesh.borderSize + 1) );
		Iterator shift = Iterator::Zeros();
		for (int k = 0; k < PdeVector::M; k++) {
			shift(s) = (dx(k) > 0) ? 1 : -1;
			for (int i = 0; i < (int)src.size(); i++) {
				src[(size_t)i] = mesh.pde(it + shift * i);
			}
			ans.setColumn(k, EqualDistanceLineInterpolator<PdeVector>::
					minMaxInterpolate(src, fabs(dx(k)) / mesh.h(s)));
		}
		return ans;
	}
	
	
};


} // namespace cubic
} // namespace gcm


#endif // LIBGCM_CUBIC_GRIDCHARACTERISTICMETHOD_HPP
