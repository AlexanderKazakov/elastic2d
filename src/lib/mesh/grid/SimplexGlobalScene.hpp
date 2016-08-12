#ifndef LIBGCM_SIMPLEXGLOBALSCENE_HPP
#define LIBGCM_SIMPLEXGLOBALSCENE_HPP

#include <lib/mesh/grid/AbstractGlobalScene.hpp>
#include <lib/mesh/grid/SimplexGrid.hpp>
#include <lib/numeric/gcm/ContactCorrector.hpp>
#include <lib/numeric/gcm/BorderCorrector.hpp>


namespace gcm {

class Engine;

/**
 * Global scene of the program -- a conglomeration of several grids --
 * for the case of calculation on simplex grids
 */
template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
class SimplexGlobalScene : public AbstractGlobalScene {
public:
	
	typedef SimplexGrid<Dimensionality, TriangulationT>    Grid;
	typedef typename Grid::Triangulation                   Triangulation;
	typedef typename Grid::Iterator                        Iterator;
	typedef typename Grid::GridId                          GridId;
	typedef typename Grid::RealD                           RealD;
	typedef typename Grid::MatrixDD                        MatrixDD;
	
	typedef typename Triangulation::VertexHandle           VertexHandle;
	typedef typename Triangulation::CellHandle             CellHandle;
	
	static const GridId EmptySpaceFlag = Grid::EmptySpaceFlag;
	
	typedef AbstractContactCorrector<Grid>                 ContactCorrector;
	typedef typename ContactCorrector::NodesContact        NodesContact;
	
	/// pair of grids in contact
	typedef std::pair<GridId, GridId>                      GridsPair;
	
	struct Contact {
		/// list of nodes pairs in contact
		std::list<NodesContact> nodesInContact;
		/// corrector to handle these contacts
		std::shared_ptr<ContactCorrector> contactCorrector;
	};
	
	
	typedef AbstractBorderCorrector<Grid>                  BorderCorrector;
	typedef typename BorderCorrector::NodeBorder           NodeBorder;
	
	struct Border {
		/// list of border nodes
		std::list<NodeBorder> borderNodes;
		/// corrector to handle these nodes
		std::shared_ptr<BorderCorrector> borderCorrector;
		/// used only at instantiation step, not in calculations
		std::shared_ptr<Area> correctionArea;
	};
	
	
	SimplexGlobalScene(const Task& task, Engine* engine_ = nullptr);
	virtual ~SimplexGlobalScene() { }
	
	
	/**
	 * Actions to perform after all grids would be constructed --
	 * find and remember all contacts and borders.
	 */
	virtual void afterGridsConstruction(const Task& task) override;
	
	
	/**
	 * Perform next time step -- stages of GCM and border/contact correction
	 */
	virtual void nextTimeStep() override;
	
	
private:
	Engine * const engine;
	
	/// global triangulation of the whole calculation space
	Triangulation triangulation;
	
	/// on/off points motion
	bool movable = false;
	
	
	/// Current basis of calculations
	MatrixDD calculationBasis = linal::randomBasis(MatrixDD());
	
	
	/// all contact conditions of all grids
	std::map<GridsPair, Contact> contacts;
	
	/// it's possible to have several border conditions for one grid
	/// all border conditions of all grids
	std::map<GridId, std::vector<Border>> borders;
	
	
	friend class SimplexGrid<Dimensionality, TriangulationT>;
	USE_AND_INIT_LOGGER("gcm.SimplexGlobalScene")
	
	
	void createNewCalculationBasis();
	void correctContactsAndBorders(const int stage);
	
	void createContacts(const Task& task);
	void createBorders(const Task& task);
	
	void addNode(const VertexHandle vh, const GridsPair gridsIds);
	void addContactNode(const VertexHandle vh, const GridsPair gridsIds);
	void addBorderNode(const VertexHandle vh, const GridId gridId);
	
	
	std::set<GridId> incidentGridsIds(const VertexHandle vh) const {
		std::list<CellHandle> incidentCells = triangulation.allIncidentCells(vh);
		std::set<GridId> ans;
		for (CellHandle ch : incidentCells) {
			ans.insert(ch->info().getGridId());
		}
		return ans;
	}
	
	
};


}

#endif // LIBGCM_SIMPLEXGLOBALSCENE_HPP
