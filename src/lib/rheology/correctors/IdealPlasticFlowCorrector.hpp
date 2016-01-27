#ifndef LIBGCM_IDEALPLASTICFLOWCORRECTOR_HPP
#define LIBGCM_IDEALPLASTICFLOWCORRECTOR_HPP

#include <cmath>
#include <lib/util/Types.hpp>

namespace gcm {
	/**
	 * The simplest plasticity flow model corrector
	 */
	class IdealPlasticFlowCorrector {
	public:
		IdealPlasticFlowCorrector(const real _yieldStrength) : yieldStrength(_yieldStrength) { };
		real yieldStrength = 0;

		template<class TNode>
		void plasticFlowCorrector(TNode& node);

	};

};

#endif // LIBGCM_IDEALPLASTICFLOWCORRECTOR_HPP
