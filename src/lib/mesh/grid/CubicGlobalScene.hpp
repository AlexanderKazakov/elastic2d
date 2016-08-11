#ifndef LIBGCM_CUBICGLOBALSCENE_HPP
#define LIBGCM_CUBICGLOBALSCENE_HPP

#include <lib/mesh/grid/AbstractGlobalScene.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>
#include <lib/Engine.hpp>


namespace gcm {


template<int Dimensionality>
struct CubicGlobalScene : public AbstractGlobalScene {
	typedef CubicGrid<Dimensionality> Grid;
	
	CubicGlobalScene(const Task&, Engine* engine_ = nullptr) :
			engine(engine_) { }
	virtual ~CubicGlobalScene() { }
	
	virtual void afterGridsConstruction(const Task&) override { }
	
	virtual void nextTimeStep() override {
		assert_true(engine);
		for (int stage = 0; stage < Dimensionality; stage++) {
			for (const auto body : engine->bodies) {
				body.second.solver->privateStage(stage, Clock::TimeStep());
			}
		}
	}
	
private:
	Engine * const engine;
	
};



}

#endif // LIBGCM_CUBICGLOBALSCENE_HPP
