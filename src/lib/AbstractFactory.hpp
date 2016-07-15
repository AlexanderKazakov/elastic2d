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
	virtual Solver* createSolver(const Task& task) const = 0;
	virtual Snapshotter* createSnapshotter(Snapshotters::T snapshotterId) const = 0;

};
template<
        template<typename, typename, typename> class TMesh,
        typename TModel, typename TGrid, typename TMaterial
        >
struct ProgramAbstactFactory : public VirtualProgramAbstactFactory {
	virtual Solver* createSolver(const Task& task) const override {
		return new DefaultSolver<TMesh<TModel, TGrid, TMaterial> >(task);
	}

	virtual Snapshotter* createSnapshotter(Snapshotters::T snapshotterId) const override {
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
	virtual VirtualProgramAbstactFactory* create() const = 0;
};
template<typename TModel, typename TGrid, typename TMaterial>
struct MeshFactory : public VirtualMeshFactory {
	virtual VirtualProgramAbstactFactory* create() const override {
		return new ProgramAbstactFactory<DefaultMesh, TModel, TGrid, TMaterial>();
	}
};


/**
 * Choose grid
 */
struct VirtualGridFactory {
	virtual VirtualMeshFactory* create(Grids::T) const = 0;
};
template<typename TModel, typename TMaterial>
struct GridFactory : public VirtualGridFactory {
	virtual VirtualMeshFactory* create(Grids::T gridId) const override {
		switch (gridId) {
			case Grids::T::CUBIC:
				return new MeshFactory<TModel, CubicGrid<TModel::DIMENSIONALITY>, TMaterial>();

		case Grids::T::CGAL:
			switch (TModel::DIMENSIONALITY) {
				case 2:
					return create<(TModel::DIMENSIONALITY == 2) &&
							(std::is_same<TMaterial, IsotropicMaterial>::value), Cgal2DGrid>();
				case 3:
					return create<(TModel::DIMENSIONALITY == 3) &&
							(std::is_same<TMaterial, IsotropicMaterial>::value), Cgal3DGrid>();
				default:
					THROW_INVALID_ARG("CGAL grids has dimensions 2 or 3 only and support only isotropic materials by now");
			}

		default:
			THROW_INVALID_ARG("The type of grid is unknown or unsuitable for specified type of model and material");
		}
	}
	
private:
	template<bool valid_use, typename GridType>
	typename std::enable_if<valid_use, VirtualMeshFactory*>::type create() const {
		return new MeshFactory<TModel, GridType, TMaterial>();
	}
	template<bool valid_use, typename GridType>
	typename std::enable_if<!valid_use, VirtualMeshFactory*>::type create() const {
		THROW_INVALID_ARG("The type of grid is not suitable for specified type of material and rheology model");
	}
};


/**
 * Choose rheology model
 */
struct VirtualModelFactory {
	virtual VirtualGridFactory* create(Models::T) const = 0;
};
template<int Dimensionality, typename TMaterial>
struct ModelFactory : public VirtualModelFactory {
	virtual VirtualGridFactory* create(Models::T modelId) const override {
		switch (modelId) {
			case Models::T::ELASTIC:
				return new GridFactory<ElasticModel<Dimensionality>, TMaterial>();
			case Models::T::ACOUSTIC:
				return new GridFactory<AcousticModel<Dimensionality>, TMaterial>();
			default:
				THROW_INVALID_ARG("Model is unknown or incompatible for specified type of material and dimensionality");
		}
	}
};
template<int Dimensionality>
struct ModelFactory<Dimensionality, OrthotropicMaterial> : public VirtualModelFactory {
	virtual VirtualGridFactory* create(Models::T modelId) const override {
		switch (modelId) {
			case Models::T::ELASTIC:
				return new GridFactory<ElasticModel<Dimensionality>, OrthotropicMaterial>();
			default:
				THROW_INVALID_ARG("Model is unknown or incompatible for specified type of material and dimensionality");
		}
	}
};


/**
 * Choose material
 */
struct VirtualMaterialFactory {
	virtual VirtualModelFactory* create(Materials::T) const = 0;
};
template<int Dimensionality>
struct MaterialFactory : public VirtualMaterialFactory {
	virtual VirtualModelFactory* create(Materials::T materialId) const override {
		switch (materialId) {
			case Materials::T::ISOTROPIC:
				return new ModelFactory<Dimensionality, IsotropicMaterial>();
			case Materials::T::ORTHOTROPIC:
				return new ModelFactory<Dimensionality, OrthotropicMaterial>();
			default:
				THROW_UNSUPPORTED("Unknown type of material");
		}
	}
};
template<>
struct MaterialFactory<1> : public VirtualMaterialFactory {
	virtual VirtualModelFactory* create(Materials::T materialId) const override {
		switch (materialId) {
			case Materials::T::ISOTROPIC:
				return new ModelFactory<1, IsotropicMaterial>();
			default:
				THROW_UNSUPPORTED("Material is unknown or incompatible for specified dimensionality");
		}
	}
};


/**
 * Choose space dimensionality
 */
struct DimensionalityFactory {
	VirtualMaterialFactory* create(const int dimensionality) const {
		switch (dimensionality) {
			case 1:
				return new MaterialFactory<1>();
			case 2:
				return new MaterialFactory<2>();
			case 3:
				return new MaterialFactory<3>();
			default:
				THROW_UNSUPPORTED("Invalid dimensionality");
		}
	}
};


/**
 * The program instantiator
 */
struct Factory {
	static VirtualProgramAbstactFactory* create(const Task& task) {
		return create(
				task.dimensionality,
				task.materialId,
				task.gridId,
				task.modelId);
	}

private:
	static VirtualProgramAbstactFactory* create(
			const int dimensionality,
			Materials::T materialId,
			Grids::T gridId,
			Models::T modelId) {
	
		return (new DimensionalityFactory())
				->create(dimensionality)
				->create(materialId)
				->create(modelId)
				->create(gridId)
				->create();
	}
};


}


#endif // LIBGCM_ABSTRACTFACTORY_HPP
