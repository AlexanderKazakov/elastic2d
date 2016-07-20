#ifndef LIBGCM_GRIDADAPTER_HPP
#define LIBGCM_GRIDADAPTER_HPP

namespace gcm {

/**
 * This is a full-featured grid type which can be used 
 * as a template parameter for meshes.
 * It shares some big basic grid with several similar adapters.
 * That allows to keep several meshes in one big geometric structure.
 * @tparam BasicGrid big common geometric structure type
 */
template<typename BasicGrid>
class GridAdapter {
public:
	
	/// Space dimensionality
	static const int DIMENSIONALITY = BasicGrid::DIMENSIONALITY;
	
	/// An *estimation* of maximal possible number of vertices connected 
	/// with some inner vertex (it can be more in a very rare cases)
	static const int MAX_NUMBER_OF_NEIGHBOR_VERTICES = 
			BasicGrid::MAX_NUMBER_OF_NEIGHBOR_VERTICES;
	
	
	/// @name Iterators 
	///@{
	typedef BasicGrid::Iterator Iterator;
	
	typedef Iterator ForwardIterator;
	ForwardIterator begin() const { return 0; }
	ForwardIterator end() const { return sizeOfRealNodes(); }
	typedef Iterator VtkIterator;
	VtkIterator vtkBegin() const { return begin(); }
	VtkIterator vtkEnd() const { return end(); }
	
	/** Iteration over all border nodes */
	///@{
	typedef typename std::vector<size_t>::const_iterator BorderIterator;
	BorderIterator borderBegin() const { return borderIndices.begin(); }
	BorderIterator borderEnd() const { return borderIndices.end(); }
	///@}
	/** Iteration over all inner nodes */
	///@{
	typedef typename std::vector<size_t>::const_iterator InnerIterator;
	InnerIterator innerBegin() const { return innerIndices.begin(); }
	InnerIterator innerEnd() const { return innerIndices.end(); }
	///@}
	
	
private:
	BasicGrid* commonGrid;
	
	GridAdapter(const Task& task);
	virtual ~GridAdapter() { }
};


}

#endif // LIBGCM_GRIDADAPTER_HPP
