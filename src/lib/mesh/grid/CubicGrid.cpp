#include <fstream>
#include <algorithm>

#include <lib/mesh/grid/CubicGrid.hpp>
#include <lib/util/task/InitialCondition.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;


CubicGrid::CubicGrid(const Task &task) :
		StructuredGrid(task),
		borderSize(task.borderSize),
		sizes(calculateSizes(task)),
		indexMaker(calculateIndexMaker(task)),
		startR(calculateStartR(task)),
		h(calculateH(task)) {

	this->minimalSpatialStep = fmin(h(0), fmin(h(1), h(2)));

	assert_ge(borderSize, 1);
	for (int j = 0; j < 3; j++) {
		if (sizes(j) != 1) assert_ge(sizes(j), borderSize);
		assert_gt(h(j), 0.0);
		assert_eq(h(j), h(j)); // this is supposed to catch NaN
		assert_eq(startR(j), startR(j));
	}
}


CubicGrid::Int3 CubicGrid::calculateSizes(const Task& task) const {
	Int3 _sizes = task.sizes;
	// MPI - we divide the grid among processes equally along x-axis
	_sizes(0) = numberOfNodesAlongXPerOneCore(task);
	// in order to keep specified in task number of nodes
	if (this->rank == this->numberOfWorkers - 1)
		_sizes(0) = task.sizes(0) - numberOfNodesAlongXPerOneCore(task) * (this->numberOfWorkers - 1);
	// MPI - (end)
	return _sizes;
}

CubicGrid::Int3 CubicGrid::calculateIndexMaker(const Task& task) const {
	Int3 _indexMaker = {0, 0, 0};
	Int3 _sizes = calculateSizes(task);
	int _borderSize = task.borderSize;

	switch (task.dimensionality) {
		case 1:
			_indexMaker(0) = 1;
			break;
		case 2:
			_indexMaker(0) = 2 * _borderSize + _sizes(1);
			_indexMaker(1) = 1;
			break;
		case 3:
			_indexMaker(0) = (2 * _borderSize + _sizes(1)) * (2 * _borderSize + _sizes(2));
			_indexMaker(1) = 2 * _borderSize + _sizes(2);
			_indexMaker(2) = 1;
			break;
		default:
			THROW_INVALID_ARG("Invalid task dimensionality");
	}

	return _indexMaker;
}

CubicGrid::Real3 CubicGrid::calculateStartR(const Task& task) const {
	Real3 _startR = task.startR;
	// MPI - we divide the grid among processes equally along x-axis
	_startR(0) += (this->rank * numberOfNodesAlongXPerOneCore(task))
				  * calculateH(task)(0); // mpi parallel along X axis
	// MPI - (end)
	return _startR;
}

CubicGrid::Real3 CubicGrid::calculateH(const Task& task) const {
	return linal::plainDivision(task.lengthes, task.sizes - Int3({1, 1, 1}));
}

int CubicGrid::numberOfNodesAlongXPerOneCore(const Task& task) const {
	return (int) std::round((real) task.sizes(0) / this->numberOfWorkers);
}


















