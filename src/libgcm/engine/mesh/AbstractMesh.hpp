#ifndef LIBGCM_ABSTRACTMESH_HPP
#define LIBGCM_ABSTRACTMESH_HPP

#include <libgcm/util/infrastructure/infrastructure.hpp>
#include <libgcm/util/Enum.hpp>

namespace gcm {

class Task;

template<typename TGrid>
class AbstractMesh : public TGrid {
public:
	typedef TGrid                           Grid;
	typedef typename Grid::GridId           GridId;
	typedef typename Grid::MatrixDD         MatrixDD;
	typedef typename Grid::ConstructionPack ConstructionPack;
	
	AbstractMesh(const GridId id_, const ConstructionPack& constructionPack) :
			Grid(id_, constructionPack) { }
	
	virtual ~AbstractMesh() { }
	
	
	virtual Models::T getModelType() const = 0;
	virtual Materials::T getMaterialType() const = 0;
	
	virtual void setUpPde(const Task& task) = 0;
	virtual real getMaximalEigenvalue() const = 0;
	
	/**
	 * (For SimplexGrid only)
	 * Change calculation basis of the inner nodes only
	 * AND change their gcm-matrices according to the basis
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


}


#endif // LIBGCM_ABSTRACTMESH_HPP
