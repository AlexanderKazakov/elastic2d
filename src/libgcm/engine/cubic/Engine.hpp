#ifndef LIBGCM_CUBIC_ENGINE_HPP
#define LIBGCM_CUBIC_ENGINE_HPP

#include <libgcm/engine/AbstractEngine.hpp>
#include <libgcm/engine/cubic/AbstractFactory.hpp>


namespace gcm {
namespace cubic {


/**
 * Engine implements grid-characteristic method 
 * with dimensional splitting (stages) on rectangular grids
 */
template<int Dimensionality>
class Engine : public AbstractEngine {
public:
	static const int DIMENSIONALITY = Dimensionality;
	typedef CubicGrid<DIMENSIONALITY>            Grid;
	typedef AbstractMesh<Grid>                   Mesh;
	typedef typename Grid::IntD                  IntD;
	typedef typename Grid::RealD                 RealD;
	typedef typename Grid::GridId                GridId;
	typedef typename Grid::ConstructionPack      GridConstructionPack;
	
	Engine(const Task& task);
	virtual ~Engine() { }
	
	
	static int numberOfNodesAlongXPerOneCore(const Task::CubicGrid& task) {
		return (int) std::round((real) task.sizes.at(0) / Mpi::Size());
	}
	
	
	/// for tests TODO replace
	const AbstractGrid* getAbstractGrid() const {
		assert_eq(bodies.size(), 1);
		return bodies.front().grid.get();
	}
	
	
protected:
	virtual void nextTimeStep() override;
	virtual real estimateTimeStep() override;
	virtual void writeSnapshots(const int step) override;
	
	
private:
	struct Body {
		std::shared_ptr<Mesh> grid;
		
		std::shared_ptr<GridCharacteristicMethodBase> gcm;
		
		typedef std::shared_ptr<AbstractOde> OdePtr;
		std::vector<OdePtr> odes;
		
		typedef std::shared_ptr<Snapshotter> SnapPtr;
		std::vector<SnapPtr> snapshotters;
	};
	
	/// list of all bodies
	std::vector<Body> bodies;
	
	
	/** 
	 * Factories of meshes, methods, odes and snapshotters
	 * for certain models and materials creation
	 */
	std::shared_ptr<AbstractFactoryBase<Grid>>
	createAbstractFactory(const Task::Body& body);
	
	
	/// Functions for GridConstructionPack creation @{
	static GridConstructionPack createGridConstructionPack(
			const Task& task, const GridId gridId);
	static IntD calculateSizes(const Task::CubicGrid& task);
	static RealD calculateH(const Task::CubicGrid& task);
	static RealD calculateStartR(const Task::CubicGrid& task);
	///@}
	
	
};


} // namespace cubic
} // namespace gcm


#endif // LIBGCM_CUBIC_ENGINE_HPP
