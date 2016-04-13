#ifndef LIBGCM_INITIALCONDITION_HPP
#define LIBGCM_INITIALCONDITION_HPP

#include "Task.hpp"

namespace gcm {

template<typename TModel, typename TMaterial>
class InitialCondition {
public:
	typedef typename TModel::PdeVariables    PdeVariables;
	typedef typename TModel::PdeVector       PdeVector;
	typedef typename TModel::GCM_MATRICES    GCM_MATRICES;

	/** Convert task terms of initial conditions to own format */
	void initialize(const Statement& statement);

	/** Apply initial conditions to node assume that its coordinates is coords */
	void apply(PdeVector& v, const Real3& coords) const;

private:
	struct PdeCondition {
		PdeCondition(std::shared_ptr<Area> area_, PdeVector pdeVector_) :
			area(area_), pdeVector(pdeVector_) { }
		std::shared_ptr<Area> area;
		PdeVector pdeVector;
	};
	std::vector<PdeCondition> pdeConditions;

};


}


#endif // LIBGCM_INITIALCONDITION_HPP
