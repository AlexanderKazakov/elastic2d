#ifndef LIBGCM_VELOCITYSIGMAVARIABLES_HPP
#define LIBGCM_VELOCITYSIGMAVARIABLES_HPP

#include <lib/linal/linal.hpp>
#include <lib/rheology/variables/GetSetter.hpp>
#include <lib/util/Enum.hpp>

namespace gcm {

/**
 * Velocity and symmetric tension tensor components
 * @tparam Dimensionality space dimensionality
 * @tparam Size size of PDE-vector (the usual value is default, but, for example,
 * \sigma_{zz} may appear in 2D - then use explicit specialization)
 */
template<int Dimensionality, 
         int Size = Dimensionality + (Dimensionality * (Dimensionality + 1) ) / 2>
struct VelocitySigmaVariables : public linal::Vector<Size> {

	typedef linal::Vector<Size>                               PdeVector;
	
	typedef GetSetter<PdeVector>                              GETSETTER;
	typedef Vector3GetSetter<PdeVector>                       VECTOR3GETSETTER;

	typedef std::map<PhysicalQuantities::T, GETSETTER>        QuantitiesMap;
	typedef std::map<PhysicalQuantities::T, VECTOR3GETSETTER> Vector3Map;

	static const int DIMENSIONALITY = Dimensionality;	
	
	using PdeVector::PdeVector;
	using PdeVector::operator=;
	
	/** Access to velocity */
	///@{
	real velocity(const int i) const {
		assert_lt(i, DIMENSIONALITY);
		return (*this)(i);
	}
	
	real& velocity(const int i) {
		assert_lt(i, DIMENSIONALITY);
		return (*this)(i);
	}
	///@}

	/** Access to sigma as symmetric matrix */
	///@{
	real sigma(const int i, const int j) const {
		assert_lt(i, DIMENSIONALITY);
		assert_lt(j, DIMENSIONALITY);
		return (*this)(DIMENSIONALITY + 
		               linal::Symmetric<DIMENSIONALITY>::getIndex(i, j));
	}

	real& sigma(const int i, const int j) {
		assert_lt(i, DIMENSIONALITY);
		assert_lt(j, DIMENSIONALITY);
		return (*this)(DIMENSIONALITY + 
		               linal::Symmetric<DIMENSIONALITY>::getIndex(i, j));
	}
	///@}

	real getPressure() const {
		real trace_ = 0;
		for (int i = 0; i < DIMENSIONALITY; i++) {
			trace_ += this->sigma(i, i);
		}
		return - trace_ / DIMENSIONALITY;
	}
	void setPressure(const real& pressure) {
		linal::clear(*this);
		for (int i = 0; i < DIMENSIONALITY; i++) {
			sigma(i, i) = - pressure;
		}
	}

	/**
	 * Second tensor deviator invariant
	 * J2 == sqrt( 0.5 * s_{ij}*s_{ij} ),
	 * where s_{ij} == sigma_{ij} + pressure * (i == j)
	 */
	real getJ2() const {
		real J22 = 0;
		real pressure = getPressure();
		for (int i = 0; i < DIMENSIONALITY; i++) {
			for (int j = 0; j < DIMENSIONALITY; j++) {
				J22 += (sigma(i, j) + (i == j) * pressure) *
					   (sigma(i, j) + (i == j) * pressure) / 2;
			}
		}
		return sqrt(J22);
	}

	/** 
	 * @name Getters and Setters
	 * @see GetSetter.hpp for explanations
	 */
	///@{
	static const QuantitiesMap QUANTITIES;
	static const Vector3Map    VECTORS;

	template<int i>
	static real GetV(const PdeVector& variablesToGetFrom) {
		static_assert(i < DIMENSIONALITY, "Index out of range");
		return static_cast<const VelocitySigmaVariables&>
				(variablesToGetFrom).velocity(i);
	}

	template<int i> 
	static void SetV(const real& value, PdeVector& variablesToSetTo) {
		static_assert(i < DIMENSIONALITY, "Index out of range");
		static_cast<VelocitySigmaVariables&>(variablesToSetTo).velocity(i) = value;
	}

	template<int i, int j>
	static real GetSigma(const PdeVector& variablesToGetFrom) {
		static_assert(i < DIMENSIONALITY && j < DIMENSIONALITY, "Index out of range");
		return static_cast<const VelocitySigmaVariables&>
				(variablesToGetFrom).sigma(i, j);
	}

	template<int i, int j>
	static void SetSigma(const real& value, PdeVector& variablesToSetTo) {
		static_assert(i < DIMENSIONALITY && j < DIMENSIONALITY, "Index out of range");
		static_cast<VelocitySigmaVariables&>(variablesToSetTo).sigma(i, j) = value;
	}

	static real GetPressure(const PdeVector& variablesToGetFrom) {
		return static_cast<const VelocitySigmaVariables&>
				(variablesToGetFrom).getPressure();
	}

	static void SetPressure(const real& value, PdeVector& variablesToSetTo) {
		static_cast<VelocitySigmaVariables&>(variablesToSetTo).setPressure(value);
	}

	static Real3 GetVelocity(const PdeVector& variablesToGetFrom) {
		Real3 ans = {0, 0, 0};
		for (int i = 0; i < DIMENSIONALITY; i++) {
			ans(i)= static_cast<const VelocitySigmaVariables&>
						(variablesToGetFrom).velocity(i);
		}
		return ans;
	}

	static void SetVelocity(const Real3& value, PdeVector& variablesToSetTo) {
		for (int i = 0; i < DIMENSIONALITY; i++) {
			static_cast<VelocitySigmaVariables&>
					(variablesToSetTo).velocity(i) = value(i);
		}
	}
	///@}
};
}


#endif // LIBGCM_VELOCITYSIGMAVARIABLES_HPP
