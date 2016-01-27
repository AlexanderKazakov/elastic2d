#include <lib/rheology/correctors/IdealPlasticFlowCorrector.hpp>

using namespace gcm;


template<class TNode>
void IdealPlasticFlowCorrector::plasticFlowCorrector(TNode &node) {

	real pressure = node.getPressure();
	real J2 = node.getJ2();
	// Correction parameter
	real x = yieldStrength / J2;

	if (x < 1) {
		for (int i = 0; i < TNode::DINENSIONALITY; i++) {
			for (int j = 0; j < TNode::DIMENSIONALITY; j++) {
				node.sigma(i, j) = x * (node.sigma(i, j) + (i == j) * pressure) - (i == j) * pressure;
			}
		}
	}
}
