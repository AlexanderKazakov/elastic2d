#ifndef LIBGCM_ENGINE_HPP
#define LIBGCM_ENGINE_HPP

#include <lib/util/task/Task.hpp>
#include <lib/mesh/grid/AbstractGlobalScene.hpp>
#include <lib/numeric/solvers/Solver.hpp>
#include <lib/util/snapshot/Snapshotter.hpp>
#include <lib/util/Logging.hpp>


/**
 * @mainpage
 *
 * GCM - grid-characteristic method
 * ================================
 * 
 * Engine for numerical simulations of wave and associated processes
 * -----------------------------------------------------------------
 * 
 * ## About
 *
 * 1-2-3-dimensional simulation of wave and associated processes
 *
 */


namespace gcm {

template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
class SimplexGlobalScene;

template<int Dimensionality>
class CubicGlobalScene;


/**
 * Main class. Responsible for the whole process of calculation
 */
class Engine {
public:
	typedef AbstractGrid::GridId GridId;
	
	
	Engine(const Task& task_);
	~Engine();
	Engine(const Engine&) = delete;
	Engine& operator=(const Engine&) = delete;
	
	
	/**
	 * Perform all calculations
	 */
	void run();
	
	
	/**
	 * Return mesh with given id
	 */
	AbstractGrid* getAbstractMesh(const GridId id) {
		return bodies.at(id).solver->getAbstractMesh();
	}
	
	
	/// for tests
	const Solver* getSolver() const {
		assert_eq(1, bodies.size());
		return bodies.begin()->second.solver;
	}
	
	
private:
	
	struct Body {
		Solver* solver;
		std::vector<Snapshotter*> snapshotters;
	};
	
	AbstractGlobalScene* globalScene;
	
	/// Bodies sorted by unique id
	std::map<GridId, Body> bodies;
	
	
	real requiredTime = 0;
	
	
	void nextTimeStep();
	
	void estimateTimeStep();
	
	
	template<int Dimensionality, 
	         template<int, typename, typename> class TriangulationT>
	friend class SimplexGlobalScene;
	
	template<int Dimensionality>
	friend class CubicGlobalScene;
	
	USE_AND_INIT_LOGGER("gcm.Engine")
};


}

#endif // LIBGCM_ENGINE_HPP
