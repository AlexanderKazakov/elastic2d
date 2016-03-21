#ifndef LIBGCM_IDEALPLASTICFLOWCORRECTOR_HPP
#define LIBGCM_IDEALPLASTICFLOWCORRECTOR_HPP

#include <cmath>
#include <lib/util/Types.hpp>
#include <lib/util/task/Task.hpp>

namespace gcm {
	/**
	 * The simplest plasticity flow model corrector
	 */
	class IdealPlasticFlowCorrector {
	public:
		static const bool NonTrivial = true;

		IdealPlasticFlowCorrector(const Statement &) { }

		template<typename TNodePtr>
		void apply(TNodePtr node) const {
			real pressure = node->pde().getPressure();
			real J2 = node->pde().getJ2();
			// Correction parameter
			real x = J2 / node->material()->yieldStrength;

			if (x > 1) {
				for (int i = 0; i < node->pde().DIMENSIONALITY; i++) {
					for (int j = 0; j < node->pde().DIMENSIONALITY; j++) {
						node->_pde().sigma(i, j) =
								(node->pde().sigma(i, j) + (i == j) * pressure) / x - (i == j) * pressure;
					}
				}
			}
		}

	};
}

#endif // LIBGCM_IDEALPLASTICFLOWCORRECTOR_HPP
