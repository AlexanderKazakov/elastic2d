#ifndef LIBGCM_ACOUSTICVARIABLES_HPP
#define LIBGCM_ACOUSTICVARIABLES_HPP

#include <libgcm/linal/linal.hpp>
#include <libgcm/rheology/variables/GetSetter.hpp>
#include <libgcm/util/Enum.hpp>

namespace gcm {

/**
 * Velocity vector components and pressure scalar
 * @tparam Dimensionality space dimensionality
 */
template<int Dimensionality>
struct AcousticVariables : public linal::Vector<Dimensionality + 1> {
	
	typedef linal::Vector<Dimensionality + 1>                 PdeVector;
	
	typedef GetSetter<PdeVector>                              GETSETTER;
	typedef Vector3GetSetter<PdeVector>                       VECTOR3GETSETTER;

	typedef std::map<PhysicalQuantities::T, GETSETTER>        QuantitiesMap;
	typedef std::map<PhysicalQuantities::T, VECTOR3GETSETTER> Vector3Map;

	static const int DIMENSIONALITY = Dimensionality;
	typedef linal::Vector<DIMENSIONALITY>                     RealD;
	
	using PdeVector::PdeVector;
	using PdeVector::operator=;
	
	/** Access to velocity */
	///@{
	RealD getVelocity() const {
		RealD ans;
		for (int i = 0; i < DIMENSIONALITY; i++) {
			ans(i) = velocity(i);
		}
		return ans;
	}
	
	void setVelocity(const RealD& orig) {
		for (int i = 0; i < DIMENSIONALITY; i++) {
			velocity(i) = orig(i);
		}
	}
	
	real velocity(const int i) const {
		assert_lt(i, DIMENSIONALITY);
		return (*this)(i);
	}
	
	real& velocity(const int i) {
		assert_lt(i, DIMENSIONALITY);
		return (*this)(i);
	}
	///@}
	
	/** Access to pressure */
	///@{
	real  pressure() const { return (*this)(DIMENSIONALITY); }
	real& pressure()       { return (*this)(DIMENSIONALITY); }
	///@}
	
	/// Shortcut to unify access to force (not velocity)
	/// part of variables with ElasticModel
	///@{
	real getSigma() const { return pressure(); }
	void setSigma(const real& orig) { pressure() = orig; }
	///@}
	
	
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
		return static_cast<const AcousticVariables&>(variablesToGetFrom).velocity(i);
	}

	template<int i> 
	static void SetV(const real& value, PdeVector& variablesToSetTo) {
		static_assert(i < DIMENSIONALITY, "Index out of range");
		static_cast<AcousticVariables&>(variablesToSetTo).velocity(i) = value;
	}

	static real GetPressure(const PdeVector& variablesToGetFrom) {
		return static_cast<const AcousticVariables&>(variablesToGetFrom).pressure();
	}

	static void SetPressure(const real& value, PdeVector& variablesToSetTo) {
		static_cast<AcousticVariables&>(variablesToSetTo).pressure() = value;
	}

	static Real3 GetVelocity(const PdeVector& variablesToGetFrom) {
		Real3 ans = {0, 0, 0};
		for (int i = 0; i < DIMENSIONALITY; i++) {
			ans(i)= static_cast<const AcousticVariables&>(variablesToGetFrom).velocity(i);
		}
		return ans;
	}

	static void SetVelocity(const Real3& value, PdeVector& variablesToSetTo) {
		for (int i = 0; i < DIMENSIONALITY; i++) {
			static_cast<AcousticVariables&>(variablesToSetTo).velocity(i) = value(i);
		}
	}
	///@}
};
}


#endif // LIBGCM_ACOUSTICVARIABLES_HPP
