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
	virtual void changeCalculationBasis(const MatrixDD& basis) = 0;
	virtual void swapPdeTimeLayers() = 0;
};


}


#endif // LIBGCM_ABSTRACTMESH_HPP
