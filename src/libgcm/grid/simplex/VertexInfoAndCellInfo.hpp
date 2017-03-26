#ifndef VERTEXINFOANDCELLINFO_HPP
#define VERTEXINFOANDCELLINFO_HPP

#include <libgcm/util/Enum.hpp>

namespace gcm {

/** Auxiliary information stored in global triangulation cells */
template<int CELL_POINTS_NUMBER>
struct CellInfoT {
	/// Indicator that no grid owns the cell (auxiliary empty cell)
	static const GridId EmptySpaceFlag = (GridId)(-1);
	
	/// global triangulation cell can belongs to the only one grid
	GridId gridId;
	void setGridId(const GridId gridId_) { gridId = gridId_; }
	GridId getGridId() const { return gridId; }
	
	/// local indices of the cell's vertices in the order
	/// the same with their pointers (VertexHandles)
	typedef size_t LocalVertexIndex;
	LocalVertexIndex localVertexIndices[CELL_POINTS_NUMBER];
};

/** Auxiliary information stored in global triangulation vertices */
typedef size_t VertexInfo;

}

#endif // VERTEXINFOANDCELLINFO_HPP
