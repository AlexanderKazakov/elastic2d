#ifndef LIBGCM_INITIALCONDITION_HPP
#define LIBGCM_INITIALCONDITION_HPP

#include "Task.hpp"

namespace gcm {

	template<class TModel>
	class InitialCondition {
	public:
		typedef typename TModel::PdeVector       PdeVector;
		typedef typename TModel::GCM_MATRICES    GCM_MATRICES;

		/** Convert task terms of initial conditions to own format */
		void initialize(const Task& task);

		/**
		 * Apply initial conditions to node assume that its coordinates is coords.
		 * We don't use smth like node.coords because some nodes don't have coords.
		 */
		void apply(PdeVector& v, const linal::Vector3& coords) const;

	private:
		struct Condition {
			Condition(std::shared_ptr<Area> _area, PdeVector _pdeVector) : 
					area(_area), pdeVector(_pdeVector) { };
			std::shared_ptr<Area> area;
			PdeVector pdeVector;
		};
		std::vector<Condition> conditions = {};

	};
}


#endif // LIBGCM_INITIALCONDITION_HPP
