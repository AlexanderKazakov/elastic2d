#ifndef ELASTIC2D_MPISOLVER_HPP
#define ELASTIC2D_MPISOLVER_HPP


#include "lib/Mesh.hpp"

class MPISolver {
public:
	MPISolver(Mesh* mesh, Mesh* newMesh);
	/**
	 * Perform calculation of the task
	 */
	void calculate();
	/**
	 * Do next stage of splitting method
	 * @param s 0 - along X-axis, 1 - along Y-axis
	 * @param timeStep time step
	 */
	void stage(const uint s, const real& timeStep);

	bool makeSnapshots = false;
	bool splittingSecondOrder = false;

private:

	Mesh* mesh;
	Mesh* newMesh;

	void exchangeNodesWithNeighbors();

};


#endif //ELASTIC2D_MPISOLVER_HPP
