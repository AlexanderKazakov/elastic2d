#include <fstream>
#include <algorithm>

#include <lib/mesh/grid/CubicGrid.hpp>
#include <lib/util/task/InitialCondition.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;


CubicGrid::CubicGrid(const Task& task) :
		StructuredGrid(task),
		borderSize(task.cubicGrid.borderSize),
		sizes(task.cubicGrid.sizes),
		indexMaker(calculateIndexMaker(task.cubicGrid)),
		startR(task.cubicGrid.startR),
		h(task.cubicGrid.h) {
	
	assert_ge(borderSize, 1);
	for (int i = 0; i < task.cubicGrid.dimensionality; i++) {
		assert_ge(sizes(i), borderSize);
		assert_gt(h(i), 0);
	}
	for (int i = task.cubicGrid.dimensionality; i < 3; i++) {
		assert_eq(sizes(i), 1);
		assert_eq(h(i), std::numeric_limits<real>::max());
		assert_eq(indexMaker(i), 0);
	}
	for (int i = 0; i < 3; i++) {
		// this is supposed to catch NaN
		assert_eq(h(i), h(i));
		assert_eq(startR(i), startR(i));
	}
	
	std::cout << Mpi::Rank() << " # " << "sizes = " << sizes << " h = " << h << "sR = " << startR;
}

void CubicGrid::preprocessTask(Task::CubicGrid& task) {
	// check before precalculations
	assert_ge(task.dimensionality, 1);
	assert_le(task.dimensionality, 3);
	assert_ge(task.borderSize, 1);
	for (int i = 0; i < task.dimensionality; i++) {
		assert_ge(task.sizes(i), task.borderSize);
		assert_ge(task.lengths(i), 0);
		assert_ge(task.h(i), 0);
		// either lengths or h can be specified
		if (task.lengths(0) > 0) { // lengths specified
			assert_eq(task.h(i), 0);
			assert_gt(task.lengths(i), 0);
		} else { // h specified
			assert_eq(task.lengths(i), 0);
			assert_gt(task.h(i), 0);
		}
	}
	for (int i = task.dimensionality; i < 3; i++) {
		task.sizes(i) = 1;
		task.lengths(i) = std::numeric_limits<real>::max();
		task.h(i) = std::numeric_limits<real>::max();
		task.startR(i) = 0;
	}

	if (task.lengths(0) > 0) { // if lengths not h specified
		task.h = linal::plainDivision(task.lengths, task.sizes - Int3({1, 1, 1}));
	}
	// with respect to MPI partition
	// TODO - remove this shit
	auto sizes_ = calculateSizes(task);
	auto startR_ = calculateStartR(task);
	task.sizes = sizes_;
	task.startR = startR_;
}

Int3 CubicGrid::calculateSizes(const Task::CubicGrid& task) {
	Int3 _sizes = task.sizes;
	if (Mpi::ForceSequence()) {
		return _sizes;
	}
	
	// MPI - we divide the grid among processes equally along x-axis
	_sizes(0) = numberOfNodesAlongXPerOneCore(task);
	// in order to keep specified in task number of nodes
	if (Mpi::Rank() == Mpi::Size() - 1) {
		_sizes(0) = task.sizes(0) - 
				numberOfNodesAlongXPerOneCore(task) * (Mpi::Size() - 1);
	}
	
	return _sizes;
}

Real3 CubicGrid::calculateStartR(const Task::CubicGrid& task) {
	Real3 _startR = task.startR;
	if (Mpi::ForceSequence()) {
		return _startR;
	}
	
	// MPI - we divide the grid among processes equally along x-axis
	_startR(0) += Mpi::Rank() * numberOfNodesAlongXPerOneCore(task) * task.h(0);
	
	return _startR;
}

int CubicGrid::numberOfNodesAlongXPerOneCore(const Task::CubicGrid& task) {
	return (int) std::round((real) task.sizes(0) / Mpi::Size());
}

Int3 CubicGrid::calculateIndexMaker(const Task::CubicGrid& task) {

	Int3 _indexMaker = {0, 0, 0};
	switch (task.dimensionality) {
		case 1:
			_indexMaker(0) = 1;
			break;
		case 2:
			_indexMaker(0) = 2 * task.borderSize + task.sizes(1);
			_indexMaker(1) = 1;
			break;
		case 3:
			_indexMaker(0) = (2 * task.borderSize + task.sizes(1)) * (2 * task.borderSize + task.sizes(2));
			_indexMaker(1) = 2 * task.borderSize + task.sizes(2);
			_indexMaker(2) = 1;
			break;
		default:
			THROW_INVALID_ARG("Invalid task dimensionality");
	}

	return _indexMaker;
}


