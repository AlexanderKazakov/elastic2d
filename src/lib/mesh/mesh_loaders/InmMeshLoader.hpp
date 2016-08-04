#ifndef LIBGCM_INMMESHLOADER_HPP
#define LIBGCM_INMMESHLOADER_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <lib/util/FileUtils.hpp>
#include <lib/util/StringUtils.hpp>
#include <lib/linal/linal.hpp>
#include <lib/util/Logging.hpp>


namespace gcm {

/**
 * Class for loading unstructured tetrahedral 3d meshes 
 * by Yuri Vassilevski group from Institute of Numerical Mathematics
 */
class InmMeshLoader {
public:
	
	static const int Dimensionality = 3;
	static const int NumberOfCellVertices = Dimensionality + 1;
	typedef linal::Vector<Dimensionality>            Point;
	typedef std::array<size_t, NumberOfCellVertices> Cell;
	typedef int                                      Material;
	
	static const char delimiter = ' ';
	static const size_t EmptyMaterialFlag = (size_t)(-1);
	
	
	/**
	 * Load triangulation from file into CGAL 3D triangulation.
	 * Loaded triangulation may slightly differ from the source one
	 * because of some special CGAL rules for triangulations
	 */
	template<typename Triangulation>
	static void load(const std::string fileName, Triangulation& triangulation) {
		
		typedef typename Triangulation::Vertex_handle  VertexHandle;
		typedef typename Triangulation::Cell_handle    CellHandle;
		
		std::vector<Real3> points;
		std::map<Cell, Material> materials;
		
		USE_AND_INIT_LOGGER("InmMeshLoader")
		LOG_INFO("Start reading from file \"" << fileName << "\" ...");
		readFromFile(fileName, points, materials);
		
		
		LOG_INFO("Start adding points ...");
		CellHandle insertHint = CellHandle();
		for (size_t i = 0; i < points.size(); i++) {
			CgalPoint point = makeCgalPoint(points[i]);
			VertexHandle vh = triangulation.insert(point, insertHint);
			vh->info() = i + 1; // number of the vertex in INM mesh format
			insertHint = vh->cell();
			if (i % 100000 == 0 && i != 0) { LOG_INFO(i << " points have been loaded"); }
		}
		
		
		LOG_INFO("Start adding materials to cells ...");
		size_t counter = 0, matchCounter = 0;
		for (auto cell  = triangulation.all_cells_begin();
		          cell != triangulation.all_cells_end(); ++cell) {
			
			cell->info().setGridId(EmptyMaterialFlag);
			if (triangulation.is_infinite(cell)) { continue; }
			
			Cell inmCell = {
					cell->vertex(0)->info(),
					cell->vertex(1)->info(),
					cell->vertex(2)->info(),
					cell->vertex(3)->info(),
			};
			std::sort(inmCell.begin(), inmCell.end());
			
			const auto material = materials.find(inmCell);
			if (material != materials.end()) {
				++matchCounter;
				cell->info().setGridId((size_t)(*material).second);
			}
			
			if (++counter % 500000 == 0) { LOG_INFO(counter << " cells have been loaded"); }
		}
		
		
		LOG_INFO("Total number of given cells: " << materials.size());
		LOG_INFO("Total number of matched cells: " << matchCounter);
		size_t missedCellsCount = materials.size() - matchCounter;
		real missedCellsRatio = (real) missedCellsCount / (real) materials.size();
		LOG_INFO("Total number of missed cells: " << missedCellsCount
				<< ", percentage: " << missedCellsRatio * 100 << "%");
		
		
		correctHangedCells(triangulation);
	}
	
	
	static void readFromFile(const std::string fileName,
			std::vector<Real3>& points,
			std::map<Cell, Material>& materials) {
		points.clear();
		materials.clear();
		
		std::ifstream input;
		FileUtils::openTextFileStream(input, fileName);
		
		readPoints(input, points);
		readCells(input, materials);
		checkEndOfFile(input);
		
		FileUtils::closeFileStream(input);
	}
	
	
private:
	
	template<typename Triangulation>
	static void correctHangedCells(Triangulation& triangulation) {
	/// change material flag for cells which neighbors 
	/// all have a different material
		
		size_t emptyHangsCounter = 0, otherHangsCounter = 0;
		USE_AND_INIT_LOGGER("InmMeshLoader")
		LOG_INFO("Start replacing hanged cells ...");
		
		for (auto cell =  triangulation.finite_cells_begin();
				  cell != triangulation.finite_cells_end(); ++cell) {
			
			size_t cellMaterialFlag = cell->info().getGridId();
			size_t neighbor0MaterialFlag = cell->neighbor(0)->info().getGridId();
			size_t neighbor1MaterialFlag = cell->neighbor(1)->info().getGridId();
			size_t neighbor2MaterialFlag = cell->neighbor(2)->info().getGridId();
			size_t neighbor3MaterialFlag = cell->neighbor(3)->info().getGridId();
			
			if ( neighbor0MaterialFlag != cellMaterialFlag &&
				 neighbor0MaterialFlag == neighbor1MaterialFlag &&
				 neighbor0MaterialFlag == neighbor2MaterialFlag &&
				 neighbor0MaterialFlag == neighbor3MaterialFlag ) {
			// if all four neighbors have the same material different from the cell one
				
				if (cellMaterialFlag == EmptyMaterialFlag) { ++emptyHangsCounter; }
				else { ++otherHangsCounter; }
				
				cell->info().setGridId(neighbor0MaterialFlag); // set material from neighbors
			}
		}
		
		LOG_INFO("Replaced " << emptyHangsCounter << " single empty cells and "
				<< otherHangsCounter << " single non-empty cells");
	}
	
	
	static void readPoints(std::ifstream& input, std::vector<Real3>& points) {
		std::string numberOfPointsStr;
		std::getline(input, numberOfPointsStr);
		assert_eq(StringUtils::split(numberOfPointsStr, delimiter).size(), 1);
		size_t numberOfPoints = std::stoul(numberOfPointsStr);
		assert_ge(numberOfPoints, NumberOfCellVertices);
		
		points.reserve(numberOfPoints);
		for (size_t i = 0; i < numberOfPoints; i++) {
			std::string pointStr;
			std::getline(input, pointStr);
			
			const auto coords = StringUtils::split(pointStr, delimiter);
			assert_eq(coords.size(), Dimensionality);
			
			real x = std::stod(coords[0]);
			real y = std::stod(coords[1]);
			real z = std::stod(coords[2]);
			points.push_back(Real3({x, y, z}));
		}
	}
	
	
	static void readCells(std::ifstream& input, std::map<Cell, Material>& materials) {
		std::string numberOfCellsStr;
		std::getline(input, numberOfCellsStr);
		assert_eq(StringUtils::split(numberOfCellsStr, delimiter).size(), 1);
		size_t numberOfCells = std::stoul(numberOfCellsStr);
		assert_ge(numberOfCells, 1);
		
		for (size_t i = 0; i < numberOfCells; i++) {
			std::string cellStr;
			std::getline(input, cellStr);
			
			const auto cellStrs = StringUtils::split(cellStr, delimiter);
			assert_eq(cellStrs.size(), NumberOfCellVertices + 1);
			
			Cell cell;
			for (size_t j = 0; j < NumberOfCellVertices; j++) {
				cell[j] = std::stoul(cellStrs[j]);
			}
			std::sort(cell.begin(), cell.end()); // sort to simplify search later
			Material material = std::stoi(cellStrs[NumberOfCellVertices]);
			
			materials.insert({cell, material});
		}
	}
	
	
	static void checkEndOfFile(std::ifstream& input) {
		std::string zeroLine; // end of file
		std::getline(input, zeroLine);
		assert_eq(StringUtils::split(zeroLine, delimiter).size(), 1);
		assert_eq(std::stoi(zeroLine), 0);
	}
	
	
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef K::Point_3 CgalPoint;
	
	static CgalPoint makeCgalPoint(const Point& p) {
		return CgalPoint(p(0), p(1), p(2));
	}
	
	
};


}


#endif // LIBGCM_INMMESHLOADER_HPP
