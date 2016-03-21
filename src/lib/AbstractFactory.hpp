#ifndef LIBGCM_ABSTRACTFACTORY_HPP
#define LIBGCM_ABSTRACTFACTORY_HPP

#include <lib/util/snapshot/snapshotters.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/rheology/models/Model.hpp>
#include <lib/util/Logging.hpp>


namespace gcm {
	
	struct AbstractFinalFactory {
		virtual Solver* createSolver(const Task& task) const = 0;
		virtual Snapshotter* createSnapshotter(Snapshotters::T snapshotterId) const = 0;
	};
	template<typename TMesh>
	struct FinalFactory : public AbstractFinalFactory {
		virtual Solver* createSolver(const Task& task) const override {
			return new DefaultSolver<TMesh>(task);
		}
		virtual Snapshotter* createSnapshotter(Snapshotters::T snapshotterId) const override {
			switch (snapshotterId) {
				case Snapshotters::T::VTK:
					return new VtkSnapshotter<TMesh>();
				case Snapshotters::T::BIN2DSEISM:
					return new Binary2DSeismograph<TMesh>();
				case Snapshotters::T::DETECTOR:
					return new Detector<TMesh>();
				default:
					THROW_UNSUPPORTED("Unknown type of snapshotter");
			}
		}
	};

	struct AbstractMeshFactory {
		virtual AbstractFinalFactory* create() const = 0;
	};
	template<typename TModel, typename TGrid, typename TMaterial>
	struct MeshFactory : public AbstractMeshFactory {
		virtual AbstractFinalFactory* create() const override {
			return new FinalFactory<DefaultMesh<TModel, TGrid, TMaterial>>();
		}
	};
	
	struct AbstractModelFactory {
		virtual AbstractMeshFactory* create(Models::T) const = 0;
	};
	template<typename TGrid, typename TMaterial>
	struct ModelFactory : public AbstractModelFactory {
		virtual AbstractMeshFactory* create(Models::T modelId) const override {
			switch (modelId) {
				case Models::T::ELASTIC1D:
					return new MeshFactory<Elastic1DModel, TGrid, TMaterial>();
				case Models::T::ELASTIC2D:
					return new MeshFactory<Elastic2DModel, TGrid, TMaterial>();
				case Models::T::ELASTIC3D:
					return new MeshFactory<Elastic3DModel, TGrid, TMaterial>();
				default:
					THROW_UNSUPPORTED("Unknown type of rheology model");
			}
		}
	};
	
	struct AbstractGridFactory {
		virtual AbstractModelFactory* create(Grids::T) const = 0;
	};
	template<typename TMaterial>
	struct GridFactory : public AbstractGridFactory {
		virtual AbstractModelFactory* create(Grids::T gridId) const override {
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
	
	struct MaterialFactory {
		AbstractGridFactory* create(Materials::T materialId) const {
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
	
	struct Factory {
		static AbstractFinalFactory* create(const Task& task) {
			return create(task.materialId, task.gridId, task.modelId);
		}
		
	private:
		static AbstractFinalFactory* create(Materials::T materialId, Grids::T gridId, Models::T modelId) {
			return (new MaterialFactory())->create(materialId)->create(gridId)->create(modelId)->create();
		}
	};
}

#endif // LIBGCM_ABSTRACTFACTORY_HPP
