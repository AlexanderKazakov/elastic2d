#ifndef LIBGCM_SIMPLEX_ENGINE_HPP
#define LIBGCM_SIMPLEX_ENGINE_HPP

#include <libgcm/engine/AbstractEngine.hpp>
#include <libgcm/grid/simplex/SimplexGrid.hpp>
#include <libgcm/engine/simplex/AbstractFactory.hpp>
#include <libgcm/engine/simplex/ContactCorrector.hpp>
#include <libgcm/engine/simplex/BorderCorrector.hpp>


namespace gcm {
namespace simplex {


/**
 * Engine implements grid-characteristic method 
 * with dimensional splitting (stages) on simplex grids
 */
template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
class Engine : public AbstractEngine {
public:
	
	typedef SimplexGrid<Dimensionality, TriangulationT>    Grid;
	typedef AbstractMesh<Grid>                             Mesh;
	typedef typename Grid::Triangulation                   Triangulation;
	typedef typename Grid::ConstructionPack                GridConstructionPack;
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
	
	
	Engine(const Task& task);
	virtual ~Engine() { }
	
	std::shared_ptr<const Mesh> getMesh(const GridId gridId) const {
		return getBody(gridId).grid;
	}
	
protected:
	
	/**
	 * Perform next time step -- stages of GCM and border/contact correction
	 */
	virtual void nextTimeStep() override;
	virtual void writeSnapshots(const int step) override;
	
	virtual real estimateTimeStep() override {
		/// minimal among all bodies
		real minimalTimeStep = std::numeric_limits<real>::max();
		for (const Body& body : bodies) {
			real bodyTimeStep = CourantNumber *
					body.grid->getAverageHeight() /
					body.grid->getMaximalEigenvalue();
			if (bodyTimeStep < minimalTimeStep) {
				minimalTimeStep = bodyTimeStep;
			}
		}
		return minimalTimeStep;
	}
	
	
private:
	struct Body {
		std::shared_ptr<Mesh> grid;
		
		std::shared_ptr<GridCharacteristicMethodBase> gcm;
		
		typedef std::shared_ptr<AbstractOde> OdePtr;
		std::vector<OdePtr> odes;
		
		typedef std::shared_ptr<Snapshotter> SnapPtr;
		std::vector<SnapPtr> snapshotters;
		
		/// it's possible to have several border conditions for one grid
		std::vector<Border> borders;
	};
	
	/// global triangulation of the whole calculation space
	Triangulation triangulation;
	
	/// list of all bodies
	std::vector<Body> bodies;
	
	/// all contact conditions of all grids
	std::map<GridsPair, Contact> contacts;
	
	/// on/off points motion
	bool movable = false;
	
	struct CalculationBasis {
	/// Current basis of calculations --
	/// i'th gcm stage is performed along i'th column of the matrix
		MatrixDD basis = MatrixDD::Zeros();
		bool createNewRandomAtEachTimeStep = false;
	} calculationBasis;
	
	friend class SimplexGrid<Dimensionality, TriangulationT>;
	USE_AND_INIT_LOGGER("gcm.simplex.Engine")
	
	
	Body& getBody(const GridId gridId) {
		for (Body& body : bodies) {
			if (body.grid->id == gridId) { return body; }
		}
		THROW_INVALID_ARG("There isn't a body with given id");
	}
	const Body& getBody(const GridId gridId) const {
		for (const Body& body : bodies) {
			if (body.grid->id == gridId) { return body; }
		}
		THROW_INVALID_ARG("There isn't a body with given id");
	}
	
	
	void initializeCalculationBasis(const Task& task) {
		std::vector<real> taskBasis = task.calculationBasis;
		calculationBasis.createNewRandomAtEachTimeStep = taskBasis.empty();
		if (calculationBasis.createNewRandomAtEachTimeStep) {
			LOG_INFO("Use new random calculation basis at each time step");
		} else {
			assert_eq(taskBasis.size(), Dimensionality * Dimensionality);
			calculationBasis.basis.copyFrom(taskBasis);
			LOG_INFO("Use constant calculation basis:" << calculationBasis.basis);
			for (const Body& body : bodies) {
				body.grid->changeCalculationBasis(calculationBasis.basis);
			}
		}
	}
	void createNewCalculationBasis() {
		if (!calculationBasis.createNewRandomAtEachTimeStep) { return; }
		calculationBasis.basis = linal::randomBasis(calculationBasis.basis);
		LOG_INFO("New calculation basis:" << calculationBasis.basis);
		for (const Body& body : bodies) {
			body.grid->changeCalculationBasis(calculationBasis.basis);
		}
	}
	
	void correctContactsAndBorders(const int stage);
	
	void createMeshes(const Task& task);
	void createContacts(const Task& task);
	
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
	
	
	/** Creation of the factory of meshes and snapshotters */
	std::shared_ptr<AbstractFactoryBase<Grid>>
	createAbstractFactory(const Task::Body& body) {
		if (body.materialId != Materials::T::ISOTROPIC) {
			THROW_UNSUPPORTED("Unsupported material type");
		}
		
		switch (body.modelId) {
			case (Models::T::ACOUSTIC):
				return std::make_shared<AbstractFactory<
						AcousticModel<Dimensionality>,
						Grid, IsotropicMaterial, DefaultMesh>>();
			case (Models::T::ELASTIC):
				return std::make_shared<AbstractFactory<
						ElasticModel<Dimensionality>,
						Grid, IsotropicMaterial, DefaultMesh>>();
			default:
				THROW_UNSUPPORTED("Unknown model type");
		}
		
	}
	
};


} // namespace simplex
} // namespace gcm


#endif // LIBGCM_SIMPLEX_ENGINE_HPP
