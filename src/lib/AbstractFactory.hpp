#ifndef LIBGCM_ABSTRACTFACTORY_HPP
#define LIBGCM_ABSTRACTFACTORY_HPP

#include <lib/util/snapshot/snapshotters.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/rheology/models/Model.hpp>
#include <lib/util/Logging.hpp>


namespace gcm {

/**
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
		return create < TModel::DIMENSIONALITY == 2 &&
		       std::is_same<TGrid, CubicGrid>::value > ();
	}

	template<bool valid_use>
	typename std::enable_if<valid_use, Snapshotter*>::type create() {
		return new Binary2DSeismograph<TMesh<TModel, TGrid, TMaterial> >();
	}

	template<bool valid_use>
	typename std::enable_if<!valid_use, Snapshotter*>::type create() {
		THROW_INVALID_ARG(
		        "The type of snapshotter is not suitable for specified type of grid and rheology model");
	}

};
template<
        template<typename, typename, typename> class TMesh,
        typename TModel, typename TGrid, typename TMaterial
        >
struct SnapshotterFactory<Detector, TMesh, TModel, TGrid, TMaterial> {
	Snapshotter* create() {
		return create<std::is_same<TGrid, CubicGrid>::value>();
	}

	template<bool valid_use>
	typename std::enable_if<valid_use, Snapshotter*>::type create() {
		return new Detector<TMesh<TModel, TGrid, TMaterial> >();
	}

	template<bool valid_use>
	typename std::enable_if<!valid_use, Snapshotter*>::type create() {
		THROW_INVALID_ARG(
		        "The type of snapshotter is not suitable for specified type of grid and rheology model");
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
			return SnapshotterFactory<VtkSnapshotter, TMesh, TModel, TGrid,
			                          TMaterial>().create();
		case Snapshotters::T::BIN2DSEISM:
			return SnapshotterFactory<Binary2DSeismograph, TMesh, TModel, TGrid,
			                          TMaterial>().create();
		case Snapshotters::T::DETECTOR:
			return SnapshotterFactory<Detector, TMesh, TModel, TGrid,
			                          TMaterial>().create();
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
 * Choose rheology model
 */
struct VirtualModelFactory {
	virtual VirtualMeshFactory* create(Models::T) const = 0;

};
template<typename TGrid, typename TMaterial>
struct ModelFactory : public VirtualModelFactory {
	virtual VirtualMeshFactory* create(Models::T modelId) const override {
		switch (modelId) {
		case Models::T::ELASTIC1D:
			return new MeshFactory<Elastic1DModel, TGrid, TMaterial>();
		case Models::T::ELASTIC2D:
			return new MeshFactory<Elastic2DModel, TGrid, TMaterial>();
		case Models::T::ELASTIC3D:
			return new MeshFactory<Elastic3DModel, TGrid, TMaterial>();
		default:
			THROW_INVALID_ARG(
			        "The type of rheology model is unknown or unsuitable for specified type of grid and material");
		}
	}

};
template<typename TGrid>
struct ModelFactory<TGrid, OrthotropicMaterial> : public VirtualModelFactory {
	virtual VirtualMeshFactory* create(Models::T modelId) const override {
		switch (modelId) {
		case Models::T::ELASTIC2D:
			return new MeshFactory<Elastic2DModel, TGrid, OrthotropicMaterial>();
		case Models::T::ELASTIC3D:
			return new MeshFactory<Elastic3DModel, TGrid, OrthotropicMaterial>();
		default:
			THROW_INVALID_ARG(
			        "The type of rheology model is unknown or unsuitable for specified type of grid and material");
		}
	}

};
template<typename TMaterial>
struct ModelFactory<Cgal2DGrid, TMaterial> : public VirtualModelFactory {
	virtual VirtualMeshFactory* create(Models::T modelId) const override {
		switch (modelId) {
		case Models::T::ELASTIC2D:
			return new MeshFactory<Elastic2DModel, Cgal2DGrid, TMaterial>();
		default:
			THROW_INVALID_ARG(
			        "The type of rheology model is unknown or unsuitable for specified type of grid and material");
		}
	}

};
template<>         // to prevent ambiguous template specialization
struct ModelFactory<Cgal2DGrid, OrthotropicMaterial> : public VirtualModelFactory {
	virtual VirtualMeshFactory* create(Models::T modelId) const override {
		switch (modelId) {
		case Models::T::ELASTIC2D:
			return new MeshFactory<Elastic2DModel, Cgal2DGrid, OrthotropicMaterial>();
		default:
			THROW_INVALID_ARG(
			        "The type of rheology model is unknown or unsuitable for specified type of grid and material");
		}
	}

};

/**
 * Choose grid
 */
struct VirtualGridFactory {
	virtual VirtualModelFactory* create(Grids::T) const = 0;

};
template<typename TMaterial>
struct GridFactory : public VirtualGridFactory {
	virtual VirtualModelFactory* create(Grids::T gridId) const override {
		switch (gridId) {
		case Grids::T::CUBIC:
			return new ModelFactory<CubicGrid, TMaterial>();
		case Grids::T::CGAL2D:
			return new ModelFactory<Cgal2DGrid, TMaterial>();
		default:
			THROW_UNSUPPORTED("Unknown type of grid");
		}
	}

};

/**
 * Choose material
 */
struct MaterialFactory {
	VirtualGridFactory* create(Materials::T materialId) const {
		switch (materialId) {
		case Materials::T::ISOTROPIC:
			return new GridFactory<IsotropicMaterial>();
		case Materials::T::ORTHOTROPIC:
			return new GridFactory<OrthotropicMaterial>();
		default:
			THROW_UNSUPPORTED("Unknown type of material");
		}
	}

};

/**
 * Create ProgramAbstractFactory to instantiate the program.
 * This is "ProgramAbstractFactory factory".
 */
struct Factory {
	static VirtualProgramAbstactFactory* create(const Task& task) {
		return create(task.materialId, task.gridId, task.modelId);
	}

	private:
		static VirtualProgramAbstactFactory* create(Materials::T materialId,
		                                            Grids::T gridId,
		                                            Models::T modelId) {
			return (new MaterialFactory())->create(materialId)->create(gridId)->create(
			               modelId)->create();
		}

};
}

#endif // LIBGCM_ABSTRACTFACTORY_HPP
