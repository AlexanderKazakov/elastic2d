#ifndef GCM_STRUCTUREDGRIDBORDERCONDITIONS_HPP
#define GCM_STRUCTUREDGRIDBORDERCONDITIONS_HPP

#include <lib/util/Concepts.hpp>
#include <lib/util/task/Task.hpp>

namespace gcm {
	template<class TNode> class StructuredGrid;

	template<class TNode>
	class StructuredGridBorderConditions {

		std::map<CUBIC_BORDERS, BorderCondition::T> borderConditions;

	public:
		void initialize(const Task& task);
		void applyBorderConditions(StructuredGrid<TNode>* mesh);

	};
}

#endif // GCM_STRUCTUREDGRIDBORDERCONDITIONS_HPP
