#ifndef LIBGCM_ABSTRACTGLOBALSCENE_HPP
#define LIBGCM_ABSTRACTGLOBALSCENE_HPP


namespace gcm {

class Task;

/**
 * Abstract base class for global scene of the program --
 * a conglomeration of several grids
 */
class AbstractGlobalScene {
public:
	virtual ~AbstractGlobalScene() { }
	
	
	/**
	 * Actions to perform after all grids would be constructed
	 */
	virtual void afterGridsConstruction(const Task&) = 0;
	
	
	/**
	 * Perform next time step (part, connected with waves propagation)
	 */
	virtual void nextTimeStep() = 0;
	
};


}

#endif // LIBGCM_ABSTRACTGLOBALSCENE_HPP
