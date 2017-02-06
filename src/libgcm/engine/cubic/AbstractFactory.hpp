#ifndef LIBGCM_CUBIC_ABSTRACTFACTORY_HPP
#define LIBGCM_CUBIC_ABSTRACTFACTORY_HPP

#include <memory>

#include <libgcm/util/snapshot/snapshotters.hpp>
#include <libgcm/rheology/ode/Ode.hpp>
#include <libgcm/engine/cubic/DefaultMesh.hpp>
#include <libgcm/engine/cubic/GridCharacteristicMethod.hpp>
#include <libgcm/engine/cubic/ContactConditions.hpp>
#include <libgcm/engine/cubic/BorderConditions.hpp>
#include <libgcm/rheology/materials/materials.hpp>


namespace gcm {
namespace cubic {

/**
 * Factory of meshes, snapshotters, numerical methods, etc
 * for the cubic gcm engine
 */
template<typename TGrid>
class AbstractFactoryBase {
public:
	typedef std::shared_ptr<AbstractMesh<TGrid>>          MeshPtr;
	typedef std::shared_ptr<GridCharacteristicMethodBase> GcmPtr;
	typedef std::shared_ptr<AbstractOde>                  OdePtr;
	typedef std::shared_ptr<Snapshotter>                  SnapPtr;
	typedef std::shared_ptr<AbstractContactCopier<TGrid>> ContactPtr;
	typedef std::shared_ptr<AbstractBorderConditions>     BorderPtr;
	
	typedef typename TGrid::ConstructionPack     GridConstructionPack;
	typedef typename TGrid::PartIterator         PartIterator;
	typedef typename AbstractGrid::GridId        GridId;
	
	
	virtual MeshPtr createMesh(const Task& task, const GridId gridId,
			const GridConstructionPack& constructionPack,
			const size_t numberOfNextPdeTimeLayers) = 0;
	
	virtual GcmPtr createGcm(const Task& task) = 0;
	
	virtual BorderPtr createBorder(const Task& task, const MeshPtr mesh) = 0;
	
	virtual OdePtr createOde(const Odes::T type) = 0;
	
	virtual SnapPtr createSnapshotter(
			const Task& task, const Snapshotters::T type) = 0;
	
	virtual ContactPtr createContact(
			const PartIterator& iterA, const PartIterator& iterB,
			const ContactConditions::T condition,
			const Models::T neighborModel,
			const Materials::T neighborMaterial) = 0;
};


template<typename TModel, typename TGrid, typename TMaterial,
         template<typename, typename, typename> class TMesh>
class AbstractFactory : public AbstractFactoryBase<TGrid> {
public:
	typedef TMesh<TModel, TGrid, TMaterial>        Mesh;
	
	typedef AbstractFactoryBase<TGrid>             Base;
	typedef typename Base::MeshPtr                 MeshPtr;
	typedef typename Base::GcmPtr                  GcmPtr;
	typedef typename Base::OdePtr                  OdePtr;
	typedef typename Base::SnapPtr                 SnapPtr;
	typedef typename Base::GridId                  GridId;
	typedef typename Base::GridConstructionPack    GridConstructionPack;
	typedef typename Base::ContactPtr              ContactPtr;
	typedef typename Base::BorderPtr               BorderPtr;
	typedef typename Base::PartIterator            PartIterator;
	
	
	virtual MeshPtr createMesh(const Task& task, const GridId gridId,
			const GridConstructionPack& constructionPack,
			const size_t numberOfNextPdeTimeLayers) override {
		return std::make_shared<Mesh>(task, gridId,
				constructionPack, numberOfNextPdeTimeLayers);
	}
	
	virtual GcmPtr createGcm(const Task& task) override {
		return std::make_shared<GridCharacteristicMethod<Mesh>>(task);
	}
	
	virtual BorderPtr createBorder(
			const Task& task, const MeshPtr mesh) override {
		return std::make_shared<BorderConditions<Mesh>>(task, *mesh);
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
			case Snapshotters::T::SLICESNAP:
				return std::make_shared<SliceSnapshotter<Mesh>>(task);
			default:
				THROW_UNSUPPORTED("Unknown or unsupported snapshotter");
		}
	}
	
	virtual ContactPtr createContact(
			const PartIterator& iterA, const PartIterator& iterB,
			const ContactConditions::T condition,
			const Models::T neighborModel,
			const Materials::T neighborMaterial) override {
		assert_true(condition == ContactConditions::T::ADHESION); // TODO
		assert_true(Mesh::ModelType == neighborModel); // TODO
		
		typedef TMesh<TModel, TGrid, IsotropicMaterial>   IsotropicMesh;
		typedef TMesh<TModel, TGrid, OrthotropicMaterial> OrthotropicMesh;
		
		switch (neighborMaterial) {
			case Materials::T::ISOTROPIC:
				return std::make_shared<ContactCopier<
						TGrid, Mesh, IsotropicMesh>>(iterA, iterB);
			case Materials::T::ORTHOTROPIC:
				return std::make_shared<ContactCopier<
						TGrid, Mesh, OrthotropicMesh>>(iterA, iterB);
			default:
				THROW_UNSUPPORTED("Unknown material");
		}
	}
	
};


} // namespace cubic
} // namespace gcm


#endif // LIBGCM_CUBIC_ABSTRACTFACTORY_HPP
