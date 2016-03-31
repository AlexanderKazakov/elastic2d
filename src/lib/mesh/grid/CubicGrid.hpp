#ifndef LIBGCM_CUBICGRID_HPP
#define LIBGCM_CUBICGRID_HPP

#include <lib/mesh/grid/StructuredGrid.hpp>
#include <lib/linal/linal.hpp>

namespace gcm {
/**
 * Non-movable structured rectangular grid.
 * Dimensionality from 1 to 3.
 * If dimensionality is 2, the grid lies in plain XY, and size by Z is 1.
 * If dimensionality is 1, the grid lies on axis X and sizes by Y and Z are 1.
 * The slowest index is always X (if going through the memory sequentially)
 */
class CubicGrid : public StructuredGrid {
public:
	struct Iterator : public Int3 {         // TODO - replace inheritance
		Int3 bounds = {0, 0, 0};
		Iterator(const Int3 start, const Int3 bounds_) {
			(*this) = start;
			bounds = bounds_;
		}

		using Int3::Matrix;
		using Int3::operator=;
		const Iterator& operator*() const { return *this; }

		protected:
			int increment(const int i) {
				(*this)(i) = ((*this)(i) + 1) % this->bounds(i);
				return (*this)(i);
			}

	};

	/** slow-X fast-Z memory access efficient iterator */
	struct ForwardIterator : public Iterator {
		using Iterator::Iterator;
		ForwardIterator& operator++() {
			!increment(2) && !increment(1) && ((*this)(0)++);
			return (*this);
		}

		ForwardIterator begin() const { return ForwardIterator({0, 0, 0}, bounds); }
		ForwardIterator end() const { return ForwardIterator({bounds(0), 0, 0}, bounds); }
	};
	ForwardIterator begin() const { return ForwardIterator({0, 0, 0}, sizes).begin(); }
	ForwardIterator end() const { return ForwardIterator({0, 0, 0}, sizes).end(); }

	/** slow-Z fast-X vtk suitable iterator */
	struct VtkIterator : public Iterator {
		using Iterator::Iterator;
		VtkIterator& operator++() {
			!increment(0) && !increment(1) && ((*this)(2)++);
			return (*this);
		}

		VtkIterator begin() const { return VtkIterator({0, 0, 0}, bounds); }
		VtkIterator end() const { return VtkIterator({0, 0, bounds(2)}, bounds); }
	};
	VtkIterator vtkBegin() const { return VtkIterator({0, 0, 0}, sizes).begin(); }
	VtkIterator vtkEnd() const { return VtkIterator({0, 0, 0}, sizes).end(); }

	/** Iteration over some rectangular part of the grid */
	struct PartIterator : public Iterator {
		ForwardIterator relativeIterator = {0, 0, 0};
		Int3 shift = {0, 0, 0};
		PartIterator(const ForwardIterator& relIter_, const Int3& shift_) :
			Iterator(relIter_ + shift_), relativeIterator(relIter_), shift(shift_) { }
		PartIterator(const Int3 start, const Int3 min, const Int3 max) :
			PartIterator(ForwardIterator(start - min, max - min), min) { }
		using Int3::operator=;
		PartIterator& operator++() {
			++relativeIterator;
			(*this) = relativeIterator + shift;
			return (*this);
		}

		PartIterator begin() const { return PartIterator(relativeIterator.begin(), shift); }
		PartIterator end() const { return PartIterator(relativeIterator.end(), shift); }
	};
	/**
	 * Iteration over slice of cube carried across specified direction through
	 * the point with specified index (index along that direction)
	 */
	PartIterator slice(const int direction, const int index) const {
		Int3 min = {0, 0, 0}; min(direction) = index;
		Int3 max = sizes; max(direction) = index + 1;
		return PartIterator(min, min, max);
	}

	/**
	 * Iteration over rectangular box of the grid
	 * from min INclusive to max EXclusive
	 */
	PartIterator box(const Int3 min, const Int3 max) const {
		return PartIterator(min, min, max);
	}

	CubicGrid(const Task& task);
	virtual ~CubicGrid() { }

	/** Read-only access to real coordinates */
	const Real3 coords(const Iterator& it) const {
		return startR + linal::plainMultiply(it, h);
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
	 */
	size_t getIndex(const Iterator& it) const {
		return getIndex(it(0), it(1), it(2));
	}

	/**
	 * @param x x index < sizes(0)
	 * @param y y index < sizes(1)
	 * @param z z index < sizes(2)
	 * @return index in std::vector
	 */
	size_t getIndex(const int x, const int y, const int z) const {
		return (size_t)
		       (indexMaker(0) * (x + borderSize) +
		        indexMaker(1) * (y + borderSize) +
		        indexMaker(2) * (z + borderSize) );
	}

	const int borderSize;  // number of virtual nodes on border
	const Int3 sizes;      // numbers of nodes along each direction
	const Int3 indexMaker; // increments of plain array index when corresponding
	// dimensional index increases by one
	const Real3 startR;    // global coordinates of the first real node of the grid
	const Real3 h;         // spatial steps along each direction

	real getMinimalSpatialStep() const {
		real minH = fmin(h(0), fmin(h(1), h(2)));
		assert_gt(minH, 0);
		return minH;
	}

	/** Checking and precalculations of the setup properties before program starts */
	static void preprocessTask(Task::CubicGrid& task);

private:
	// functions for setup precalculations
	static Int3 calculateSizes(const Task::CubicGrid& task);
	static Real3 calculateStartR(const Task::CubicGrid& task);
	static int numberOfNodesAlongXPerOneCore(const Task::CubicGrid& task);
	static Int3 calculateIndexMaker(const Task::CubicGrid& task);

	USE_AND_INIT_LOGGER("gcm.CubicGrid")
};


}

#endif // LIBGCM_CUBICGRID_HPP
