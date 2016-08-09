#ifndef LIBGCM_CUBICGRID_HPP
#define LIBGCM_CUBICGRID_HPP

#include <lib/mesh/grid/AbstractGlobalScene.hpp>
#include <lib/mesh/grid/StructuredGrid.hpp>
#include <lib/linal/linal.hpp>
#include <lib/GlobalVariables.hpp>

namespace gcm {

class Engine;


/**
 * Non-movable structured rectangular grid.
 * Dimensionality from 1 to 3.
 * If dimensionality is 2, the grid lies in plane XY.
 * If dimensionality is 1, the grid lies on axis X.
 * The slowest index is always X (if going through the memory sequentially)
 * @tparam Dimensionality 1D, 2D or 3D
 */
template<int Dimensionality>
class CubicGrid : public StructuredGrid {
public:
	
	static const int DIMENSIONALITY = Dimensionality;
	
	typedef linal::Vector<DIMENSIONALITY>                         RealD;
	typedef linal::VectorInt<DIMENSIONALITY>                      IntD;
	
	typedef IntD                                                  Iterator;
	typedef linal::SlowXFastZ<DIMENSIONALITY>                     ForwardIterator;
	typedef linal::SlowZFastX<DIMENSIONALITY>                     VtkIterator;
	typedef linal::BoxIterator<DIMENSIONALITY, linal::SlowXFastZ> PartIterator;
	
	/// Unique number of the grid among other grids of this type
	typedef size_t GridId;
	
	struct GlobalScene : public AbstractGlobalScene {
		GlobalScene(const Task&, Engine* = nullptr) { }
		virtual ~GlobalScene() { }
		virtual void afterGridsConstruction(const Task&) override { }
		virtual void correctContacts() override { }
	};
	
	/** @name memory efficient (sequential) iteration */
	/// @{
	ForwardIterator begin() const { return ForwardIterator::begin(sizes); }
	ForwardIterator end() const { return ForwardIterator::end(sizes); }
	/// @}
	
	/** @name iteration suitable for vtk snapshotters */
	/// @{
	VtkIterator vtkBegin() const { return VtkIterator::begin(sizes); }
	VtkIterator vtkEnd() const { return VtkIterator::end(sizes); }
	/// @}
	
	/**
	 * Iteration over slice of cube carried across specified direction through
	 * the point with specified index (index along that direction)
	 */
	PartIterator slice(const int direction, const int index) const {
		assert_lt(direction, DIMENSIONALITY);
		assert_lt(index, sizes(direction));
		IntD min = IntD::Zeros(); min(direction) = index;
		IntD max = sizes; max(direction) = index + 1;
		return PartIterator(min, min, max);
	}
	
	/**
	 * Iteration over rectangular box of the grid
	 * from min INclusive to max EXclusive
	 */
	PartIterator box(const IntD min, const IntD max) const {
		return PartIterator(min, min, max);
	}
	
	CubicGrid(const Task& task, GlobalScene* globalScene, const GridId gridId_);
	virtual ~CubicGrid() { }
	
	/** Read-only access to real coordinates */
	const RealD coordsD(const Iterator& it) const {
		return startR + linal::plainMultiply(it, h);
	}
	
	/** Read-only access to real coordinates */
	const Real3 coords(const Iterator& it) const {
		Real3 ans = Real3::Zeros();
		RealD coordsD_ = coordsD(it);
		for (int i = 0; i < DIMENSIONALITY; i++) {
			ans(i) = coordsD_(i);
		}
		return ans;
	}
	
	size_t sizeOfRealNodes() const {
		return (size_t) linal::directProduct(sizes);
	}
	
	size_t sizeOfAllNodes() const {
		return (size_t) (indexMaker(0) * (2 * borderSize + sizes(0)));
	}
	
	/**
	 * @param it begin() <= iterator < end()
	 * @return index in std::vector
	 * @note - no linal::dotProduct here for performance
	 */
	size_t getIndex(const Iterator& it) const {
		size_t ans = 0;
		for (int i = 0; i < DIMENSIONALITY; i++) {
			ans += (size_t) indexMaker(i) * (size_t)(it(i) + borderSize);
		}
		return ans;
	}
	
	
	
	const int borderSize;  ///< number of virtual nodes on border
	const IntD sizes;      ///< numbers of nodes along each direction
	const IntD indexMaker; ///< increments of plain array index when corresponding
	                       ///< dimensional index increases by one
	const RealD h;         ///< spatial steps along each direction
	const RealD startR;    ///< global coordinates of the first real node of the grid
	
	GlobalScene * const globalScene;
	
	real getMinimalSpatialStep() const {
		real ans = std::numeric_limits<real>::max();
		for (int i = 0; i < DIMENSIONALITY; i++) {
			if (ans > h(i)) { ans = h(i); }
		}		
		assert_gt(ans, 0);
		return ans;
	}

	static int numberOfNodesAlongXPerOneCore(const Task::CubicGrid& task) {
		return (int) std::round((real) task.sizes.at(0) / Mpi::Size());
	}

private:
	/** functions for constructor */
	///@{ 
	IntD calculateSizes(const Task::CubicGrid& task) const;
	IntD calculateIndexMaker() const;
	RealD calculateH(const Task::CubicGrid& task) const;
	RealD calculateStartR(const Task::CubicGrid& task) const;
	///@}

	USE_AND_INIT_LOGGER("gcm.CubicGrid")
};


template<int Dimensionality>
CubicGrid<Dimensionality>::
CubicGrid(const Task& task, GlobalScene* gS, const GridId gridId_) :
		StructuredGrid(task, Grids::T::CUBIC, DIMENSIONALITY, gridId_),
		borderSize(task.cubicGrid.borderSize),
		sizes(calculateSizes(task.cubicGrid)),
		indexMaker(calculateIndexMaker()),
		h(calculateH(task.cubicGrid)),
		startR(calculateStartR(task.cubicGrid)),
		globalScene(gS) {
// note: the order in the initialization list above is important!
	
	assert_gt(borderSize, 0);
	for (int i = 0; i < DIMENSIONALITY; i++) {
		assert_ge(sizes(i), borderSize);
		assert_gt(h(i), 0);
	}
}


template<int Dimensionality>
typename CubicGrid<Dimensionality>::IntD CubicGrid<Dimensionality>::
calculateSizes(const Task::CubicGrid& task) const {
	
	IntD _sizes;
	if (!task.h.empty() && !task.lengths.empty()) {
		assert_true(task.sizes.empty());
		assert_eq(task.h.size(), DIMENSIONALITY);
		assert_eq(task.lengths.size(), DIMENSIONALITY);
		
		RealD taskLength;
		taskLength.copyFrom(task.lengths);
		RealD taskH;
		taskH.copyFrom(task.h);
	
		_sizes = linal::plainDivision(taskLength, taskH) + IntD::Ones();
	
	} else {
		assert_eq(task.sizes.size(), DIMENSIONALITY);
		_sizes.copyFrom(task.sizes);
	}
	
	if (Mpi::ForceSequence() || task.forceSequence) {
		return _sizes;
	}

	// MPI - we divide the grid among processes equally along x-axis
	_sizes(0) = numberOfNodesAlongXPerOneCore(task);
	if (Mpi::Rank() == Mpi::Size() - 1) {
	// in order to keep specified in task number of nodes
		_sizes(0) = task.sizes.at(0) -
		            numberOfNodesAlongXPerOneCore(task) * (Mpi::Size() - 1);
	}

	return _sizes;
}


template<int Dimensionality>
typename CubicGrid<Dimensionality>::IntD CubicGrid<Dimensionality>::
calculateIndexMaker() const {

	IntD _indexMaker = IntD::Zeros();
	switch (DIMENSIONALITY) {
		case 1:
			_indexMaker(0) = 1;
			break;
		case 2:
			_indexMaker(0) = 2 * borderSize + sizes(1);
			_indexMaker(1) = 1;
			break;
		case 3:
			_indexMaker(0) = (2 * borderSize + sizes(1)) * (2 * borderSize + sizes(2));
			_indexMaker(1) =  2 * borderSize + sizes(2);
			_indexMaker(2) =  1;
			break;
		default:
			THROW_INVALID_ARG("Invalid DIMENSIONALITY");
	}

	return _indexMaker;
}


template<int Dimensionality>
typename CubicGrid<Dimensionality>::RealD CubicGrid<Dimensionality>::
calculateH(const Task::CubicGrid& task) const {

	RealD _h;
	if (!task.lengths.empty() && !task.sizes.empty()) {
		assert_true(task.h.empty());
		assert_eq(task.sizes.size(), DIMENSIONALITY);
		assert_eq(task.lengths.size(), DIMENSIONALITY);
		
		RealD taskLength; 
		taskLength.copyFrom(task.lengths);
		IntD taskSizes; 
		taskSizes.copyFrom(task.sizes);
		
		_h = linal::plainDivision(taskLength, taskSizes - IntD::Ones());
	
	} else {
		assert_eq(task.h.size(), DIMENSIONALITY);
		_h.copyFrom(task.h);
	}
	
	return _h;
}


template<int Dimensionality>
typename CubicGrid<Dimensionality>::RealD CubicGrid<Dimensionality>::
calculateStartR(const Task::CubicGrid &task) const {
	RealD _startR = RealD::Zeros();
	
	if (!task.startR.empty()) {
		assert_eq(task.startR.size(), DIMENSIONALITY);
		_startR.copyFrom(task.startR);
	}
	
	if (Mpi::ForceSequence() || task.forceSequence) {
		return _startR;
	}

	// MPI - divide the grid among processes equally along x-axis
	_startR(0) += Mpi::Rank() * numberOfNodesAlongXPerOneCore(task) * h(0);

	return _startR;
}


}

#endif // LIBGCM_CUBICGRID_HPP
