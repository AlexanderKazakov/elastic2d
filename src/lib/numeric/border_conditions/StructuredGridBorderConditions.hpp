#ifndef GCM_STRUCTUREDGRIDBORDERCONDITIONS_HPP
#define GCM_STRUCTUREDGRIDBORDERCONDITIONS_HPP

#include <lib/util/Concepts.hpp>
#include <lib/util/task/Task.hpp>

namespace gcm {
	template<class TGrid>
	class StructuredGridBorderConditions {
		std::map<CUBIC_BORDERS, BorderCondition::T> borderConditions;

	public:
		void initialize(const Task& task);
		void applyBorderConditions(TGrid* mesh);

	};
}

#endif // GCM_STRUCTUREDGRIDBORDERCONDITIONS_HPP
