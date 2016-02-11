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

		void initialize(const Task &task) {
			yieldStrength = task.yieldStrength;
			assert_gt(yieldStrength, 0.0);
		};

		real yieldStrength = 0;

		template<typename TNode>
		void apply(TNode &node) const {
			real pressure = node.getPressure();
			real J2 = node.getJ2();
			// Correction parameter
			real x = J2 / yieldStrength;

			if (x > 1) {
				for (int i = 0; i < TNode::DIMENSIONALITY; i++) {
					for (int j = 0; j < TNode::DIMENSIONALITY; j++) {
						node.sigma(i, j) = (node.sigma(i, j) + (i == j) * pressure) / x - (i == j) * pressure;
					}
				}
			}
		};

	};
};

#endif // LIBGCM_IDEALPLASTICFLOWCORRECTOR_HPP
