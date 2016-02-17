#include <lib/grid/Cgal2DGrid.hpp>

using namespace gcm;


void Cgal2DGrid::initializeImpl(const Task &task) {
	LOG_INFO("Start initialization");
	mesh.initialize(task);
	minimalSpatialStep = mesh.getMinimalSpatialStep();

	initializeImplImpl(task);
}
