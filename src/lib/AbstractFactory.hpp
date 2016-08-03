#ifndef LIBGCM_ABSTRACTFACTORY_HPP
#define LIBGCM_ABSTRACTFACTORY_HPP

#include <lib/util/snapshot/snapshotters.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/rheology/models/models.hpp>
#include <lib/util/Logging.hpp>


namespace gcm {

/**
 * @file
 * The approach to instantiate the meshes, solvers, snapshotters etc,
 * templated by variety of rheology models, grids and materials etc.
 * The reason for partial template specializations is that not all existing
 * types are suitable to each other, and we need compile-time differentiation.
 * The best way to read this file is from the bottom to up.
 */

/**
 * Create snapshotter
 */
template<
        template<typename> class TSnapshotter,
        template<typename, typename, typename> class TMesh,
        typename TModel, typename TGrid, typename TMaterial
        >
struct SnapshotterFactory {
	Snapshotter* create() {
		return new TSnapshotter<TMesh<TModel, TGrid, TMaterial> >();
	}

};
template<
        template<typename, typename, typename> class TMesh,
        typename TModel, typename TGrid, typename TMaterial
        >
struct SnapshotterFactory<Binary2DSeismograph, TMesh, TModel, TGrid, TMaterial> {
	Snapshotter* create() {
		return create<std::is_same<TGrid, CubicGrid<2>>::value >();
	}

	template<bool valid_use>
	typename std::enable_if<valid_use, Snapshotter*>::type create() {
		return new Binary2DSeismograph<TMesh<TModel, TGrid, TMaterial> >();
	}
	
	template<bool valid_use>
	typename std::enable_if<!valid_use, Snapshotter*>::type create() {
		THROW_INVALID_ARG("The type of snapshotter is not suitable for specified type of grid and rheology model");
	}

};
template<
        template<typename, typename, typename> class TMesh,
        typename TModel, typename TGrid, typename TMaterial
        >
struct SnapshotterFactory<Detector, TMesh, TModel, TGrid, TMaterial> {
	Snapshotter* create() {
		return create<std::is_same<TGrid, CubicGrid<1>>::value ||
		              std::is_same<TGrid, CubicGrid<2>>::value ||
		              std::is_same<TGrid, CubicGrid<3>>::value>();
	}

	template<bool valid_use>
	typename std::enable_if<valid_use, Snapshotter*>::type create() {
		return new Detector<TMesh<TModel, TGrid, TMaterial> >();
	}

	template<bool valid_use>
	typename std::enable_if<!valid_use, Snapshotter*>::type create() {
		THROW_INVALID_ARG("The type of snapshotter is not suitable for specified type of grid and rheology model");
	}

};
template<
        template<typename, typename, typename> class TMesh,
        typename TModel, typename TGrid, typename TMaterial
        >
struct SnapshotterFactory<SliceSnapshotter, TMesh, TModel, TGrid, TMaterial> {
	Snapshotter* create() {
		return create<std::is_same<TGrid, CubicGrid<3>>::value &&
		              std::is_same<TModel, ElasticModel<3>>::value>();
	}

	template<bool valid_use>
	typename std::enable_if<valid_use, Snapshotter*>::type create() {
		return new SliceSnapshotter<TMesh<TModel, TGrid, TMaterial> >();
	}

	template<bool valid_use>
	typename std::enable_if<!valid_use, Snapshotter*>::type create() {
		THROW_INVALID_ARG("The type of snapshotter is not suitable for specified type of grid and rheology model");
	}

};


/**
 * Abstract factory to instantiate modules of the program.
 * Create solvers, snapshotters, etc.
 */
struct VirtualProgramAbstactFactory {
	virtual Solver* createSolver(const Task&, const Task::Body&,
			AbstractGlobalScene*) = 0;
	virtual Snapshotter* createSnapshotter(Snapshotters::T snapshotterId) = 0;
	
};
template<
        template<typename, typename, typename> class TMesh,
        typename TModel, typename TGrid, typename TMaterial
        >
struct ProgramAbstactFactory : public VirtualProgramAbstactFactory {
	virtual Solver* createSolver(const Task& task, const Task::Body& body,
			AbstractGlobalScene* globalScene) override {
		return new DefaultSolver<TMesh<TModel, TGrid, TMaterial>>(
				task, globalScene, body.id);
	}
	
	virtual Snapshotter* createSnapshotter(Snapshotters::T snapshotterId) override {
		switch (snapshotterId) {
			case Snapshotters::T::VTK:
				return SnapshotterFactory<VtkSnapshotter, TMesh, TModel, TGrid, TMaterial>().create();
			case Snapshotters::T::BIN2DSEISM:
				return SnapshotterFactory<Binary2DSeismograph, TMesh, TModel, TGrid, TMaterial>().create();
			case Snapshotters::T::DETECTOR:
				return SnapshotterFactory<Detector, TMesh, TModel, TGrid, TMaterial>().create();
			case Snapshotters::T::SLICESNAP:
				return SnapshotterFactory<SliceSnapshotter, TMesh, TModel, TGrid, TMaterial>().create();
			default:
				THROW_UNSUPPORTED("Unknown type of snapshotter");
		}
	}
	
};


/**
 * Choose type of mesh
 */
struct VirtualMeshFactory {
	virtual VirtualProgramAbstactFactory* create() = 0;
};
template<typename TModel, typename TGrid, typename TMaterial>
struct MeshFactory : public VirtualMeshFactory {
	virtual VirtualProgramAbstactFactory* create() override {
		return new ProgramAbstactFactory<DefaultMesh, TModel, TGrid, TMaterial>();
	}
};


/**
 * Choose rheology model
 */
struct VirtualModelFactory {
	virtual VirtualMeshFactory* create(Models::T) = 0;
};
template<typename TGrid, typename TMaterial>
struct ModelFactory : public VirtualModelFactory {
	static const int D = TGrid::DIMENSIONALITY;
	virtual VirtualMeshFactory* create(Models::T modelId) override {
		switch (modelId) {
			case Models::T::ELASTIC:
				return new MeshFactory<ElasticModel<D>, TGrid, TMaterial>();
			case Models::T::ACOUSTIC:
				return new MeshFactory<AcousticModel<D>, TGrid, TMaterial>();
			default:
				THROW_UNSUPPORTED("Unknown or incompatible type of model");
		}
	}
};
template<typename TGrid>
struct ModelFactory<TGrid, OrthotropicMaterial> : public VirtualModelFactory {
	static const int D = TGrid::DIMENSIONALITY;
	virtual VirtualMeshFactory* create(Models::T modelId) override {
		switch (modelId) {
			case Models::T::ELASTIC:
				return new MeshFactory<ElasticModel<D>, TGrid, OrthotropicMaterial>();
			default:
				THROW_UNSUPPORTED("Unknown or incompatible type of model");
		}
	}
};


/**
 * Choose type of material
 */
struct VirtualMaterialFactory {
	virtual VirtualModelFactory* create(Materials::T) = 0;
};
template<typename TGrid>
struct MaterialFactory : public VirtualMaterialFactory {
	virtual VirtualModelFactory* create(Materials::T materialId) override {
		switch (materialId) {
			case Materials::T::ISOTROPIC:
				return new ModelFactory<TGrid, IsotropicMaterial>();
			case Materials::T::ORTHOTROPIC:
				return new ModelFactory<TGrid, OrthotropicMaterial>();
			default:
				THROW_UNSUPPORTED("Unknown or incompatible type of material");
		}
	}
};
template<int D, template<int, typename, typename> class Triangulation>
struct MaterialFactory<SimplexGrid<D, Triangulation>> : public VirtualMaterialFactory {
	typedef SimplexGrid<D, Triangulation> Grid;
	virtual VirtualModelFactory* create(Materials::T materialId) override {
		switch (materialId) {
			case Materials::T::ISOTROPIC:
				return new ModelFactory<Grid, IsotropicMaterial>();
			default:
				THROW_UNSUPPORTED("Unknown or incompatible type of material");
		}
	}
};
template<>
struct MaterialFactory<CubicGrid<1>> : public VirtualMaterialFactory {
	virtual VirtualModelFactory* create(Materials::T materialId) override {
		switch (materialId) {
			case Materials::T::ISOTROPIC:
				return new ModelFactory<CubicGrid<1>, IsotropicMaterial>();
			default:
				THROW_UNSUPPORTED("Unknown or incompatible type of material");
		}
	}
};


/**
 * Choose type of grid
 */
struct VirtualGridFactory {
	virtual VirtualMaterialFactory* create(Grids::T) = 0;
	virtual AbstractGlobalScene* createGlobalScene(const Task& task) = 0;
};
template<int Dimensionality>
struct GridFactory : public VirtualGridFactory {
	virtual VirtualMaterialFactory* create(Grids::T gridId) override {
		switch (gridId) {
			case Grids::T::CUBIC:
				return new MaterialFactory<CubicGrid<Dimensionality>>();
			case Grids::T::SIMPLEX:
				return new MaterialFactory<SimplexGrid<Dimensionality, CgalTriangulation>>();
			default:
				THROW_UNSUPPORTED("Unknown or incompatible type of grid");
		}
	}
	virtual AbstractGlobalScene* createGlobalScene(const Task& task) override {
		switch (task.globalSettings.gridId) {
			case Grids::T::CUBIC:
				return new typename CubicGrid<Dimensionality>::GlobalScene(task);
			case Grids::T::SIMPLEX:
				return new typename SimplexGrid<Dimensionality, CgalTriangulation>::GlobalScene(task);
			default:
				THROW_UNSUPPORTED("Unknown or incompatible type of grid");
		}
	}
};
template<>
struct GridFactory<1> : public VirtualGridFactory {
	virtual VirtualMaterialFactory* create(Grids::T gridId) override {
		switch (gridId) {
			case Grids::T::CUBIC:
				return new MaterialFactory<CubicGrid<1>>();
			default:
				THROW_UNSUPPORTED("Unknown or incompatible type of grid");
		}
	}
	virtual AbstractGlobalScene* createGlobalScene(const Task& task) override {
		switch (task.globalSettings.gridId) {
			case Grids::T::CUBIC:
				return new typename CubicGrid<1>::GlobalScene(task);
			default:
				THROW_UNSUPPORTED("Unknown or incompatible type of grid");
		}
	}
};


/**
 * Choose space dimensionality
 */
struct DimensionalityFactory {
	VirtualGridFactory* create(const int dimensionality) {
		switch (dimensionality) {
			case 1:
				return new GridFactory<1>();
			case 2:
				return new GridFactory<2>();
			case 3:
				return new GridFactory<3>();
			default:
				THROW_UNSUPPORTED("Invalid dimensionality");
		}
	}
};


/**
 * The program instantiator
 */
struct Factory {
	
	static VirtualProgramAbstactFactory* create(
			const Task& task, const Task::Body& body) {
		return (new DimensionalityFactory())
				->create(task.globalSettings.dimensionality)
				->create(task.globalSettings.gridId)
				->create(body.materialId)
				->create(body.modelId)
				->create();
	}
	
	static AbstractGlobalScene* createGlobalScene(const Task& task) {
		return (new DimensionalityFactory())
				->create(task.globalSettings.dimensionality)
				->createGlobalScene(task);
	}
	
};


}


#endif // LIBGCM_ABSTRACTFACTORY_HPP
