#ifndef LIBGCM_INITIALCONDITION_HPP
#define LIBGCM_INITIALCONDITION_HPP

#include "Task.hpp"

namespace gcm {

	template<class TModel>
	class InitialCondition {
	public:
		typedef typename TModel::Node Node;
		typedef typename TModel::Node::Vector Vector;
		typedef typename TModel::GCM_MATRICES GCM_MATRICES;

		/** Convert task terms of initial conditions to own format */
		void initialize(const Task& task, std::shared_ptr<GCM_MATRICES> gcmMatrices);

		/**
		 * Apply initial conditions to node assume that its coordinates is coords.
		 * We don't use smth like node.coords because nodes of structured grids don't have that.
		 */
		void apply(Node& node, const linal::Vector3& coords) const;

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
