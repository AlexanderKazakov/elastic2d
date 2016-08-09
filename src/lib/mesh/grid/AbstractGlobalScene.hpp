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
	 * Apply contact correctors
	 */
	virtual void correctContacts() = 0;
};


}

#endif // LIBGCM_ABSTRACTGLOBALSCENE_HPP
