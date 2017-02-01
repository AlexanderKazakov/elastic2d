#include <gtest/gtest.h>

#include <libgcm/engine/simplex/Engine.hpp>
#include <libgcm/grid/simplex/cgal/CgalTriangulation.hpp>
#include <libgcm/util/math/Area.hpp>
#include <libgcm/rheology/models/models.hpp>

#include <libgcm/engine/mesh/DefaultMesh.hpp>
#include <libgcm/grid/simplex/SimplexGrid.hpp>


using namespace gcm;
using namespace gcm::simplex;


struct Wrapper {
	typedef DefaultMesh<AcousticModel<2>, SimplexGrid<2, CgalTriangulation>, IsotropicMaterial> Mesh;
	typedef Engine<2, CgalTriangulation> ENGINE;
	static std::shared_ptr<const Mesh> getMesh(
			const ENGINE& engine, const size_t id = 0) {
		auto grid = engine.getMesh(id);
		auto mesh = std::dynamic_pointer_cast<const Mesh>(grid);
		assert_true(mesh);
		return mesh;
	}
};


TEST(Engine, ZeroInitialization) {
	Task task;
	task.globalSettings.dimensionality = 2;
	task.globalSettings.gridId = Grids::T::SIMPLEX;
	task.globalSettings.CourantNumber = 1;
	task.globalSettings.numberOfSnaps = 10;
	task.globalSettings.stepsPerSnap = 1;
	
	task.bodies = {{1, {Materials::T::ISOTROPIC, Models::T::ACOUSTIC, {}}}};
	task.simplexGrid.spatialStep = 1.15;
	Task::SimplexGrid::Body::Border bodyBorder = {{0, 3}, {4, 0}, {0, 0}};
	task.simplexGrid.bodies = {Task::SimplexGrid::Body({1, bodyBorder, {} })};
	
	task.materialConditions.type = Task::MaterialCondition::Type::BY_BODIES;
	const auto material = std::make_shared<IsotropicMaterial>(4, 2, 1, 0, 0, 0, 0);
	task.materialConditions.byBodies.bodyMaterialMap = { {1, material} };
	
	Wrapper::ENGINE engine(task);
	engine.run();
	auto mesh = Wrapper::getMesh(engine, 1);
	for (auto it = mesh->begin(); it != mesh->end(); ++it) {
		ASSERT_EQ(mesh->pde(it), mesh->pde(it).Zeros());
		ASSERT_EQ(mesh->pdeNew(0, it), mesh->pde(it).Zeros());
		ASSERT_EQ(mesh->pdeNew(1, it), mesh->pde(it).Zeros());
	}
	
	
	Task::BorderCondition borderConditionAll;
	borderConditionAll.area = std::make_shared<InfiniteArea>();
	borderConditionAll.type = BorderConditions::T::FIXED_FORCE;
	borderConditionAll.values = {[] (real) { return 0; }};
	task.borderConditions = {borderConditionAll};
	
	Wrapper::ENGINE engine2(task);
	engine2.run();
	mesh = Wrapper::getMesh(engine2, 1);
	for (auto it = mesh->begin(); it != mesh->end(); ++it) {
		ASSERT_EQ(mesh->pde(it), mesh->pde(it).Zeros());
		ASSERT_EQ(mesh->pdeNew(0, it), mesh->pde(it).Zeros());
		ASSERT_EQ(mesh->pdeNew(1, it), mesh->pde(it).Zeros());
	}
}


