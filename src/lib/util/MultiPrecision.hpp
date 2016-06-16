#ifndef LIBGCM_MULTIPRECISION_HPP
#define LIBGCM_MULTIPRECISION_HPP

#include <cmath>
#include <mpreal.h>

#include <lib/linal/linal.hpp>

namespace gcm {

/**
 * Interface to mpfr::mpreal class. 
 * mpfr::mpreal is a cpp wrapper around gnu mpfr library.
 * @link http://www.holoborodko.com/pavel/mpfr/
 * In order to install just place header from the site to path.
 */
class MultiPrecision : public mpfr::mpreal {
public:
	typedef mpfr::mpreal Base;
	using Base::Base;
	using Base::operator=;
	
	MultiPrecision() : Base() { }
	MultiPrecision(const Base& base) : Base(base) { }
	
	/** Multiprecision matrix */
	template<int TM, int TN> using Matrix = linal::MATRIX<TM, TN, MultiPrecision>;
	/** Multiprecision vector */
	template<int TM> using Vector = linal::VECTOR<TM, MultiPrecision>;
	
	/** 
	 * Create a MultiPrecision number with given DECIMAL precision.
	 * @note be careful that the constructor MultiPrecision(T, long) 
	 * expect NOT DECIMAL but BINARY precision!
	 */
	template<typename T>
	static MultiPrecision create(const T& t, const size_t decimalPrecision) {
		return MultiPrecision(t, mpfr::digits2bits((int)decimalPrecision));
	}
	
	MultiPrecision& setDecimalPrecision(const int decimalPrecision) {
		this->setPrecision((int)mpfr::digits2bits(decimalPrecision));
		return (*this);
	}
	
	static void setDefaultDecimalPrecision(const int decimalPrecision) {
		Base::set_default_prec(mpfr::digits2bits(decimalPrecision));	
	}
	
	
	
//	operator real() const {
//#ifdef LIBGCM_DOUBLE_PRECISION
//		return toDouble();
//#else
//		return toFloat();
//#endif	
//	}
	
};


} // namespace gcm


namespace std {
	/// make gcm::MultiPrecision to be floating point and arithmetic cpp-type
	/// i.e std::is_arithmetic<gcm::MultiPrecision>::value == true
	template<> struct __is_floating_point_helper<gcm::MultiPrecision> :
			public true_type { };
} // namespace std




#endif /* LIBGCM_MULTIPRECISION_HPP */
