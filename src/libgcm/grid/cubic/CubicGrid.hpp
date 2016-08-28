#ifndef LIBGCM_CUBICGRID_HPP
#define LIBGCM_CUBICGRID_HPP

#include <libgcm/linal/linal.hpp>
#include <libgcm/engine/GlobalVariables.hpp>
#include <libgcm/grid/AbstractGrid.hpp>


namespace gcm {

/**
 * Non-movable structured rectangular grid.
 * Dimensionality from 1 to 3.
 * If dimensionality is 2, the grid lies in plane XY.
 * If dimensionality is 1, the grid lies on axis X.
 * The slowest index is always X (if going through the memory sequentially)
 * @tparam Dimensionality 1D, 2D or 3D
 */
template<int Dimensionality>
class CubicGrid : public AbstractGrid {
public:
	
	static const int DIMENSIONALITY = Dimensionality;
	
	typedef linal::Matrix<DIMENSIONALITY, DIMENSIONALITY>         MatrixDD;
	typedef linal::Vector<DIMENSIONALITY>                         RealD;
	typedef linal::VectorInt<DIMENSIONALITY>                      IntD;
	
	typedef IntD                                                  Iterator;
	typedef linal::SlowXFastZ<DIMENSIONALITY>                     ForwardIterator;
	typedef linal::SlowZFastX<DIMENSIONALITY>                     VtkIterator;
	typedef linal::BoxIterator<DIMENSIONALITY, linal::SlowXFastZ> PartIterator;
	
	
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
	
	
	/** Struct for grid constructor */
	struct ConstructionPack {
		int borderSize;
		IntD sizes;
		RealD h;
		RealD startR;
	};
	
	CubicGrid(const GridId id_, const ConstructionPack& constructionPack);
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
	
	
	real getMinimalSpatialStep() const {
		real ans = std::numeric_limits<real>::max();
		for (int i = 0; i < DIMENSIONALITY; i++) {
			if (ans > h(i)) { ans = h(i); }
		}		
		assert_gt(ans, 0);
		return ans;
	}
	
	
	
private:
	IntD calculateIndexMaker() const;
	
};



template<int Dimensionality>
CubicGrid<Dimensionality>::
CubicGrid(const GridId id_, const ConstructionPack& constructionPack) :
		AbstractGrid(id_),
		borderSize(constructionPack.borderSize),
		sizes(constructionPack.sizes),
		indexMaker(calculateIndexMaker()),
		h(constructionPack.h),
		startR(constructionPack.startR) {
// note: the order in the initialization list above is important!
	
	assert_gt(borderSize, 0);
	for (int i = 0; i < DIMENSIONALITY; i++) {
		assert_ge(sizes(i), borderSize);
		assert_gt(h(i), 0);
	}
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


}

#endif // LIBGCM_CUBICGRID_HPP
