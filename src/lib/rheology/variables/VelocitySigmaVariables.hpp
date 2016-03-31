#ifndef LIBGCM_VELOCITYSIGMAVARIABLES_HPP
#define LIBGCM_VELOCITYSIGMAVARIABLES_HPP

#include <lib/linal/linal.hpp>
#include <lib/rheology/variables/GetSetter.hpp>
#include <lib/util/Enum.hpp>

namespace gcm {

/**
 * Velocity and symmetric tension tensor components
 * @tparam Dimensionality space dimensionality
 */
template<int Dimensionality>
struct VelocitySigmaVariables {
	static const int DIMENSIONALITY = Dimensionality;
	static const int SIZE = DIMENSIONALITY + (DIMENSIONALITY * (DIMENSIONALITY + 1) ) / 2;
	union {
		real values[SIZE];
		struct {
			/** Velocity */
			real V[DIMENSIONALITY];
			/** Symmetric tension tensor components */
			real S[(DIMENSIONALITY * (DIMENSIONALITY + 1) ) / 2];
		};
	};

	/** Access to sigma as symmetric matrix */
	real sigma(const int i, const int j) const {
		return S[linal::SymmetricMatrix < Dimensionality > ::getIndex(i, j)];
	}

	real& sigma(const int i, const int j) {
		return S[linal::SymmetricMatrix < Dimensionality > ::getIndex(i, j)];
	}

	real getPressure() const;
	void setPressure(const real& pressure);

	/**
	 * Second tensor deviator invariant
	 * J2 == sqrt( 0.5 * s_{ij}*s_{ij} ),
	 * where s_{ij} == sigma_{ij} + pressure * (i == j)
	 */
	real getJ2() const;

	/**
	 * Look at GetSetter.hpp for explanations about the code below
	 */
	static const std::map<PhysicalQuantities::T,
	                      GetSetter<VelocitySigmaVariables<Dimensionality> > > QUANTITIES;
	static const std::map<PhysicalQuantities::T,
	                      Vector3GetSetter<VelocitySigmaVariables<Dimensionality> > > VECTORS;

	template<int i> static real GetV(
	        const VelocitySigmaVariables<Dimensionality>& variablesToGetFrom) {
		static_assert(i < Dimensionality, "Index out of range");
		return variablesToGetFrom.V[i];
	}

	template<int i> static void SetV(const real& value,
	                                 VelocitySigmaVariables<Dimensionality>& variablesToSetTo)
	{
		static_assert(i < Dimensionality, "Index out of range");
		variablesToSetTo.V[i] = value;
	}

	template<int i, int j> static real GetSigma(
	        const VelocitySigmaVariables<Dimensionality>& variablesToGetFrom) {
		static_assert(i < Dimensionality && j < Dimensionality, "Index out of range");
		return variablesToGetFrom.sigma(i, j);
	}

	template<int i, int j> static void SetSigma(
	        const real& value, VelocitySigmaVariables<Dimensionality>& variablesToSetTo) {
		static_assert(i < Dimensionality && j < Dimensionality, "Index out of range");
		variablesToSetTo.sigma(i, j) = value;
	}

	static real GetPressure(const VelocitySigmaVariables<Dimensionality>& variablesToGetFrom) {
		return variablesToGetFrom.getPressure();
	}

	static void SetPressure(const real& value,
	                        VelocitySigmaVariables<Dimensionality>& variablesToSetTo) {
		variablesToSetTo.setPressure(value);
	}

	static Real3 GetVelocity(const VelocitySigmaVariables<Dimensionality>& variablesToGetFrom)
	{
		Real3 ans = {0, 0, 0};
		for (int i = 0; i < Dimensionality; i++) {
			ans(i) = variablesToGetFrom.V[i];
		}
		return ans;
	}

	static void SetVelocity(const Real3& value,
	                        VelocitySigmaVariables<Dimensionality>& variablesToSetTo) {
		for (int i = 0; i < Dimensionality; i++) {
			variablesToSetTo.V[i] = value(i);
		}
	}

};
}


#endif // LIBGCM_VELOCITYSIGMAVARIABLES_HPP
