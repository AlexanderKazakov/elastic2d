#ifndef LIBGCM_DUMMYCORRECTOR_HPP
#define LIBGCM_DUMMYCORRECTOR_HPP

#include <lib/util/task/Task.hpp>

namespace gcm {
	/**
	 * Dummy (stub) corrector for the case when there aren't correctors
	 */
	class DummyCorrector {
	public:
		static const bool NonTrivial = false;

		void initialize(const Task &task) {
			SUPPRESS_WUNUSED(task);
		};

		template<typename TNode>
		void apply(TNode& node) const {
			SUPPRESS_WUNUSED(node);
		};
	};
};

#endif // LIBGCM_DUMMYCORRECTOR_HPP
