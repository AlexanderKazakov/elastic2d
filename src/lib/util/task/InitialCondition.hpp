#ifndef LIBGCM_INITIALCONDITION_HPP
#define LIBGCM_INITIALCONDITION_HPP

#include "Task.hpp"

namespace gcm {

	template<class TNode>
	class InitialCondition {
	public:
		typedef typename TNode::Vector Vector;
		typedef typename TNode::GcmMatrices GcmMatrices;

		/** Convert task terms of initial conditions to own format */
		void initialize(const Task& task);

		/**
		 * Apply initial conditions to node assume that its coordinates is coords.
		 * We don't use smth like node.coords because some nodes don't have coords.
		 */
		void apply(TNode& node, const linal::Vector3& coords) const;

	private:
		struct Condition {
			Condition(std::shared_ptr<Area> area, Vector vector) : area(area), vector(vector) { };
			std::shared_ptr<Area> area;
			Vector vector;
		};
		std::vector<Condition> conditions = {};

	};
}


#endif // LIBGCM_INITIALCONDITION_HPP
