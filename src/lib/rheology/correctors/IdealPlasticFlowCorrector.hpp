#ifndef LIBGCM_IDEALPLASTICFLOWCORRECTOR_HPP
#define LIBGCM_IDEALPLASTICFLOWCORRECTOR_HPP

#include <cmath>
#include <lib/util/Types.hpp>
#include <lib/util/task/Task.hpp>

namespace gcm {
	/**
	 * The simplest plasticity flow model corrector
	 */
	template<class TGrid>
	class IdealPlasticFlowCorrector {
	public:
		void initialize(const Task &task) {
			yieldStrength = task.material.yieldStrength;
			enable = task.plasticityFlowCorrector;
		};

		void apply(TGrid *grid) const {
			if (enable) {
				assert_gt(yieldStrength, 0);
				for (auto &node : grid->nodes)
					apply(node);
			}
		};

	protected:
		bool enable = false;
		real yieldStrength = 0;

		void apply(typename TGrid::NODE &node) const {
			real pressure = node.u.getPressure();
			real J2 = node.u.getJ2();
			// Correction parameter
			real x = J2 / yieldStrength;

			if (x > 1) {
				for (int i = 0; i < TGrid::DIMENSIONALITY; i++) {
					for (int j = 0; j < TGrid::DIMENSIONALITY; j++) {
						node.u.sigma(i, j) = (node.u.sigma(i, j) + (i == j) * pressure) / x - (i == j) * pressure;
					}
				}
			};

		};

	};
};

#endif // LIBGCM_IDEALPLASTICFLOWCORRECTOR_HPP
