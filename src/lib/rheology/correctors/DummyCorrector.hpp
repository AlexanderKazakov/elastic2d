#ifndef LIBGCM_DUMMYCORRECTOR_HPP
#define LIBGCM_DUMMYCORRECTOR_HPP

#include <lib/util/task/Task.hpp>

namespace gcm {
	/**
	 * Dummy (stub) node state corrector.
	 */
	class DummyCorrector {
	public:
		static const bool NonTrivial = false;

		void beforeStatement(const Statement&) { }

		/**
		 * This function should't be called at all, this is just stub
		 */
		template<typename TNodePtr>
		void apply(TNodePtr) const { }
	};
}

#endif // LIBGCM_DUMMYCORRECTOR_HPP
