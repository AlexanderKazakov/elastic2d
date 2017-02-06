#ifndef LIBGCM_SIMPLEX_ABSTRACTFACTORY_HPP
#define LIBGCM_SIMPLEX_ABSTRACTFACTORY_HPP

#include <memory>

#include <libgcm/util/snapshot/snapshotters.hpp>
#include <libgcm/rheology/ode/Ode.hpp>
#include <libgcm/engine/simplex/AbstractMesh.hpp>
#include <libgcm/engine/simplex/GridCharacteristicMethod.hpp>
#include <libgcm/engine/simplex/DefaultMesh.hpp>


namespace gcm {
namespace simplex {

/**
 * Factory of meshes, snapshotters, numerical methods, etc
 * for the simplex gcm engine
 */
template<typename TGrid>
class AbstractFactoryBase {
public:
	typedef std::shared_ptr<AbstractMesh<TGrid>>          MeshPtr;
	typedef std::shared_ptr<GridCharacteristicMethodBase> GcmPtr;
	typedef std::shared_ptr<AbstractOde>                  OdePtr;
	typedef std::shared_ptr<Snapshotter>                  SnapPtr;
	
	typedef typename TGrid::ConstructionPack     GridConstructionPack;
	typedef typename AbstractGrid::GridId        GridId;
	
	
	virtual MeshPtr createMesh(const Task& task, const GridId gridId,
			const GridConstructionPack& constructionPack,
			const size_t numberOfNextPdeTimeLayers) = 0;
	
	virtual GcmPtr createGcm() = 0;
	
	virtual OdePtr createOde(const Odes::T type) = 0;
	
	virtual SnapPtr createSnapshotter(const Task& task, const Snapshotters::T type) = 0;
	
};


template<typename TModel, typename TGrid, typename TMaterial>
class AbstractFactory : public AbstractFactoryBase<TGrid> {
public:
	typedef gcm::simplex::DefaultMesh<TModel, TGrid, TMaterial> Mesh;
	
	typedef AbstractFactoryBase<TGrid>             Base;
	typedef typename Base::MeshPtr                 MeshPtr;
	typedef typename Base::GcmPtr                  GcmPtr;
	typedef typename Base::OdePtr                  OdePtr;
	typedef typename Base::SnapPtr                 SnapPtr;
	typedef typename Base::GridId                  GridId;
	typedef typename Base::GridConstructionPack    GridConstructionPack;
	
	
	virtual MeshPtr createMesh(const Task& task, const GridId gridId,
			const GridConstructionPack& constructionPack,
			const size_t numberOfNextPdeTimeLayers) override {
		return std::make_shared<Mesh>(
				task, gridId, constructionPack, numberOfNextPdeTimeLayers);
	}
	
	virtual GcmPtr createGcm() override {
		return std::make_shared<GridCharacteristicMethod<Mesh>>();
	}
	
	virtual OdePtr createOde(const Odes::T type) override {
		assert_true(Odes::T::MAXWELL_VISCOSITY == type); // TODO
		return std::make_shared<MaxwellViscosityOde<Mesh>>();
	}
	
	virtual SnapPtr createSnapshotter(
			const Task& task, const Snapshotters::T type) override {
		switch (type) {
			case Snapshotters::T::VTK:
				return std::make_shared<VtkSnapshotter<Mesh>>(task);
			default:
				THROW_UNSUPPORTED("Unknown or unsupported snapshotter");
		}
	}
	
};


} // namespace simplex
} // namespace gcm


#endif // LIBGCM_SIMPLEX_ABSTRACTFACTORY_HPP
