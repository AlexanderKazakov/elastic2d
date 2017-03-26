#ifndef LIBGCM_CUBIC_ABSTRACTMESH_HPP
#define LIBGCM_CUBIC_ABSTRACTMESH_HPP

#include <libgcm/util/infrastructure/infrastructure.hpp>
#include <libgcm/util/Enum.hpp>

namespace gcm {
class Task;

namespace cubic {

/**
 * Interface from cubic::Engine to meshes of different models/materials
 * @tparam TGrid instantiation of gcm::CubicGrid
 */
template<typename TGrid>
class AbstractMesh : public TGrid {
public:
	typedef TGrid                           Grid;
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
	 * It is not called in constructor
	 */
	virtual void setUpPde(const Task& task) = 0;
	
	/** Maximal in absolute value */
	virtual real getMaximalEigenvalue() const = 0;
	
	/**
	 * Swap current PDE time layer (which is single) and next PDE time layer
	 * (which may not be single) by index 'indexOfNextPde'
	 * Useful for directional splitting into PRODUCT of stages,
	 * i.e u_{n+1} = A_1 * A_2 * A_3 * u_{n}
	 */
	virtual void swapCurrAndNextPdeTimeLayer(const int indexOfNextPde) = 0;
};

} // namespace cubic
} // namespace gcm

#endif // LIBGCM_CUBIC_ABSTRACTMESH_HPP
