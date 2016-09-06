#ifndef LIBGCM_CUBIC_ENGINE_HPP
#define LIBGCM_CUBIC_ENGINE_HPP

#include <libgcm/engine/AbstractEngine.hpp>
#include <libgcm/engine/cubic/AbstractFactory.hpp>
#include <libgcm/engine/cubic/ContactConditions.hpp>
#include <libgcm/engine/cubic/BorderConditions.hpp>


namespace gcm {
namespace cubic {


/**
 * Engine implements grid-characteristic method 
 * with dimensional splitting (stages) on rectangular grids
 */
template<int Dimensionality>
class Engine : public AbstractEngine {
public:
	static const int DIMENSIONALITY = Dimensionality;
	typedef CubicGrid<DIMENSIONALITY>            Grid;
	typedef AbstractMesh<Grid>                   Mesh;
	typedef typename Grid::IntD                  IntD;
	typedef typename Grid::RealD                 RealD;
	typedef typename Grid::AABB                  AABB;
	typedef typename Grid::GridId                GridId;
	typedef typename Grid::ConstructionPack      GridConstructionPack;
	
	Engine(const Task& task);
	virtual ~Engine() { }
	
	
	std::shared_ptr<const Mesh> getMesh(const GridId gridId) const {
		return getBody(gridId).mesh;
	}
	
	
protected:
	virtual void nextTimeStep() override;
	virtual real estimateTimeStep() override;
	virtual void writeSnapshots(const int step) override;
	
	
private:
	struct Body {
		/// For creation of all entities listed below
		std::shared_ptr<AbstractFactoryBase<Grid>> factory;
		
		std::shared_ptr<Mesh> mesh;
		
		std::shared_ptr<GridCharacteristicMethodBase> gcm;
		
		std::shared_ptr<AbstractBorderConditions> border;
		
		struct Contact {
			GridId neighborId;
			int direction;
			std::shared_ptr<AbstractContactCopier<Grid>> copier;
		};
		std::vector<Contact> contacts;
		
		typedef std::shared_ptr<AbstractOde> OdePtr;
		std::vector<OdePtr> odes;
		
		typedef std::shared_ptr<Snapshotter> SnapPtr;
		std::vector<SnapPtr> snapshotters;
		
		bool operator==(const Body& other) const {
			return mesh->id == other.mesh->id;
		}
		bool operator!=(const Body& other) const {
			return !(*this == other);
		}
	};
	
	/// list of all bodies
	std::vector<Body> bodies;
	
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
	
	
	/** 
	 * Factories of meshes, methods, odes and snapshotters
	 * for certain models and materials creation
	 */
	std::shared_ptr<AbstractFactoryBase<Grid>>
	createAbstractFactory(const Task::Body& body);
	
	
	void createGridsAndContacts(const Task& task);
	
	
	static GridConstructionPack createGridConstructionPack(
			const Task::CubicGrid& task, const GridId gridId);
	
	
};


} // namespace cubic
} // namespace gcm


#endif // LIBGCM_CUBIC_ENGINE_HPP
