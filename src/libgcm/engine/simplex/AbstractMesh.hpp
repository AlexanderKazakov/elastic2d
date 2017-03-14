#ifndef LIBGCM_SIMPLEX_ABSTRACTMESH_HPP
#define LIBGCM_SIMPLEX_ABSTRACTMESH_HPP

#include <libgcm/util/infrastructure/infrastructure.hpp>
#include <libgcm/util/Enum.hpp>

namespace gcm {

class Task;

namespace simplex {

template<typename TGrid>
class AbstractMesh : public TGrid {
public:
	typedef TGrid                           Grid;
	typedef typename Grid::GridId           GridId;
	typedef typename Grid::RealD            RealD;
	typedef typename Grid::MatrixDD         MatrixDD;
	typedef typename Grid::ConstructionPack ConstructionPack;
	
	/**
	 * Constructor: grid creation only. We delay with PDE setup
	 * (and thus setUpPde() must be called explicitly before any calculations)
	 * because it can be useful for MPI,
	 * when meshes are not completely initialized on every core
	 */
	AbstractMesh(const GridId id_, const ConstructionPack& constructionPack) :
			Grid(id_, constructionPack) { }
	virtual ~AbstractMesh() { }
	
	virtual Models::T getModelType() const = 0;
	virtual Materials::T getMaterialType() const = 0;
	
	/**
	 * Setting up the part connected with PDE: storages, matrices, etc.
	 * It is not called in constructor!
	 * @param innerBasis calculation basis for all inner nodes
	 * @param borderCalcMode affects what gcm-matrices will be on borders
	 */
	virtual void setUpPde(const Task& task, const MatrixDD& innerBasis,
			const BorderCalcMode borderCalcMode) = 0;
	
	/** Maximal in absolute value */
	virtual real getMaximalEigenvalue() const = 0;
	
	/**
	 * Change calculation basis of the inner nodes only
	 * AND change their gcm-matrices according to the basis.
	 * Border and contact nodes stay away because their basis is local.
	 * @note by now, we suppose mesh homogenity!
	 */
	virtual void setInnerCalculationBasis(const MatrixDD& basis) = 0;
	
	/**
	 * Swap current PDE time layer (which is single) and next PDE time layer
	 * (which is not single) by index 'indexOfNextPde'
	 * Useful for directional splitting into PRODUCT of stages,
	 * i.e u_{n+1} = A_1 * A_2 * A_3 * u_{n}
	 */
	virtual void swapCurrAndNextPdeTimeLayer(const int indexOfNextPde) = 0;
	
	/**
	 * Average all next PDE time layers into current PDE time layer
	 * (which is single)
	 * Useful for directional splitting into SUMM of stages,
	 * i.e u_{n+1} = (A_1 * u_{n} + A_2 * u_{n} + A_3 * u_{n}) / 3
	 */
	virtual void averageNewPdeLayersToCurrent() = 0;
};

} // namespace simplex 
} // namespace gcm

#endif // LIBGCM_SIMPLEX_ABSTRACTMESH_HPP
