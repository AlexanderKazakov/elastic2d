#ifndef ELASTIC2D_MESH_HPP
#define ELASTIC2D_MESH_HPP

#include <memory>

#include "lib/Node.hpp"
#include "lib/Interpolator.hpp"
#include "lib/Task.hpp"

class Mesh {
private:

	/* ------------------ Properties and conditions ------------------ */
	
	uint accuracyOrder = 0; // order of accuracy of spatial interpolation

	uint X = 0; // number of nodes along x direction
	uint Y = 0; // number of nodes along y direction

	real h[2] = { 0.0, /* x spatial step */
	              0.0  /* y spatial step */ };

	real tau = 0.0; // time step
	real T = 0.0; // required time

	InitialConditions initialConditions = InitialConditions::Zero;
	
	/* ------------------ Properties and conditions (end) ------------------ */

	/* PDEMatrices that is common for majority of nodes */
	std::shared_ptr<PDEMatrices> defaultMatrix;

	/* Data storage. Real nodes with assistant nodes on borders */
	Node* nodes = nullptr;

	/**
	 * Operator to iterate relatively to real nodes.
	 * Read / write access
	 * @param y y index < %Y
	 * @param x x index < %X
	 */
	inline Node& operator()(const uint y, const uint x) {
		return nodes[(2 * accuracyOrder + X) * (y + accuracyOrder) + (x + accuracyOrder)];
	};
	/**
	 * Operator to iterate relatively real nodes.
	 * Read only access
	 * @param y y index < %Y
	 * @param x x index < %X
	 */
	inline const Node& get(const uint y, const uint x) const {
		return nodes[(2 * accuracyOrder + X) * (y + accuracyOrder) + (x + accuracyOrder)];
	};


	/* Spatial interpolator */
	Interpolator interpolator;

public:
	/* Mesh factory */
	void initialize(const Task& task);

	/**
	 * Interpolate nodal values in specified points.
	 * Interpolated value for k-th point in vector %dx are
	 * stored in k-th column of returned Matrix.
	 * @param stage 0 - along X direction, 1 - along Y direction
	 * @param y y-index of the reference node
	 * @param x x-index of the reference node
	 * @param dx Vector of distances from reference node on which
	 * values should be interpolated
	 * @return Matrix with interpolated nodal values in columns
	 */
	Matrix interpolateValuesAround(const uint stage, const uint y, const uint x,
	                               const Vector& dx) const;
	/* Place in %src nodal values which are necessary for
	 * interpolation in specified point. The number of placed
	 * in values is equal to src.size()
	 */
	void findSourcesForInterpolation(const uint stage, const uint y, const uint x,
	                                 const real& dx, std::vector<Vector>& src) const;

	/* Write data to vtk file */
	void snapshot(uint step) const;


friend class SequenceSolver;

	/* ---------------- For testing purposes ---------------- */
public:
	const real& getTauForTest() const { return tau; };
	const real& getH0ForTest() const { return h[0]; };
	const real& getH1ForTest() const { return h[1]; };
	const real& getTForTest() const { return T; };

	const Node& getNodeForTest(const uint y, const uint x) const { return get(y, x); };
	const std::shared_ptr<const PDEMatrices> getDefaultMatrixForTest() const { return defaultMatrix; };

	/* Change rheology in some area
	 *
	 * @param rho2rho0 = (rho in the area) / (default rho)
	 * @param lambda2lambda0 = (lambda in the area) / (default lambda)
	 * @param mu2mu0 = (mu in the area) / (default mu)
	 */
	void changeRheology(const real& rho2rho0, const real& lambda2lambda0, const real& mu2mu0);

	/* ---------------- For testing purposes (end) ---------------- */

private:
	void applyInitialConditions();
};


#endif //ELASTIC2D_MESH_HPP
