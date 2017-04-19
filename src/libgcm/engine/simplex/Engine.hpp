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
		bool useForMulticontactNodes;
	};
	
	
	Engine(const Task& task);
	virtual ~Engine() { }
	
	std::shared_ptr<const Mesh> getMesh(const GridId gridId) const {
		return getBody(gridId).mesh;
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
//			real h = (body.mesh->getAverageHeight() +
//			          body.mesh->getMinimalHeight()) / 2;
			real h = body.mesh->getAverageHeight();
			real bodyTimeStep = CourantNumber * h /
					body.mesh->getMaximalEigenvalue();
			if (bodyTimeStep < minimalTimeStep) {
				minimalTimeStep = bodyTimeStep;
			}
		}
		return minimalTimeStep;
	}
	
	
private:
	struct Body {
		std::shared_ptr<Mesh> mesh;
		
		std::shared_ptr<GridCharacteristicMethodBase> gcm;
		
		typedef std::shared_ptr<AbstractOde> OdePtr;
		std::vector<OdePtr> odes;
		
		typedef std::shared_ptr<Snapshotter> SnapPtr;
		std::vector<SnapPtr> snapshotters;
		
		/// it's possible to have several border conditions for one mesh
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
	
	/// method of border/contacts calculation
	const BorderCalcMode borderCalcMode;
	
	/// Type of gcm-method to use for calculations
	const GcmType gcmType;
	
	/// Type of splitting by directions approach @{
	const SplittingType splittingType;
	
	/// Map between stage number and index of next PDE time layer for that stage calculation
	const std::array<int, Dimensionality> stageVsLayerMap;
	
	static std::array<int, Dimensionality> createStageVsLayerMap(
			const SplittingType splittingType) {
		std::array<int, Dimensionality> ans;
		switch(splittingType) {
			case SplittingType::PRODUCT: /// all stages to 0'th
				for (size_t i = 0; i < Dimensionality; i++) {
					ans[i] = 0;
				}
				break;
			case SplittingType::SUMM: /// all stages to each own
				for (size_t i = 0; i < Dimensionality; i++) {
					ans[i] = int(i);
				}
				break;
			default:
				THROW_BAD_CONFIG("Unknown splitting type");
		}
		return ans;
	}
	
	size_t numberOfNextPdeTimeLayers() const {
		return (size_t) stageVsLayerMap[Dimensionality - 1] + 1;
	}
	/// @}
	
	
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
			if (body.mesh->id == gridId) { return body; }
		}
		THROW_INVALID_ARG("There isn't a body with given id");
	}
	const Body& getBody(const GridId gridId) const {
		for (const Body& body : bodies) {
			if (body.mesh->id == gridId) { return body; }
		}
		THROW_INVALID_ARG("There isn't a body with given id");
	}
	
	
	void initializeCalculationBasis(const Task& task) {
		calculationBasis.createNewRandomAtEachTimeStep =
				task.calculationBasis.empty();
		if (calculationBasis.createNewRandomAtEachTimeStep) {
			calculationBasis.basis = linal::randomBasis(calculationBasis.basis);
			LOG_INFO("Use new random calculation basis at each time step");
		} else {
			assert_eq(task.calculationBasis.size(), Dimensionality * Dimensionality);
			calculationBasis.basis.copyFrom(task.calculationBasis);
			LOG_INFO("Use constant calculation basis:" << calculationBasis.basis);
		}
	}
	
	void changeCalculationBasis() {
		if (!calculationBasis.createNewRandomAtEachTimeStep) { return; }
		calculationBasis.basis = linal::randomBasis(calculationBasis.basis);
//		LOG_INFO("New calculation basis:" << calculationBasis.basis);
		for (const Body& body : bodies) {
			body.mesh->setInnerCalculationBasis(calculationBasis.basis);
		}
	}
	
	void gcmStage(const int stage, const real currentTime, const real timeStep);
	void correctContactsAndBorders(const int stage, const real timeAtNextLayer);
	void applyPlainBorderContactCorrection(const real timeForBorderCondition);
	
	void createMeshes(const Task& task);
	void createContacts(const Task& task);
	
	void addBorderOrContact(const VertexHandle vh);
	void addContactNode(const VertexHandle vh, const GridsPair gridsIds);
	void addBorderNode(const VertexHandle vh, const GridId gridId);
	
	
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
						Grid, IsotropicMaterial>>();
			case (Models::T::ELASTIC):
				return std::make_shared<AbstractFactory<
						ElasticModel<Dimensionality>,
						Grid, IsotropicMaterial>>();
			default:
				THROW_UNSUPPORTED("Unknown model type");
		}
		
	}
	
};


} // namespace simplex
} // namespace gcm


#endif // LIBGCM_SIMPLEX_ENGINE_HPP
