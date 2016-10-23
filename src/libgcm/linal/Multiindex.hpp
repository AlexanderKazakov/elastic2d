#ifndef LIBGCM_LINAL_MULTIINDEX_HPP
#define LIBGCM_LINAL_MULTIINDEX_HPP

#include <libgcm/linal/Matrix.hpp>
#include <libgcm/linal/functions.hpp>

namespace gcm {
namespace linal {

/// @name Multiindex iterators
/// @{


/** Base for multiindex iterators */
template<int D>
struct Multiindex : public VectorInt<D> {
	typedef VectorInt<D> IntD;
	
	Multiindex(const IntD& start, const IntD& bounds_) :
			IntD(start), bounds(bounds_) { }
	
	const Multiindex& operator*() const { return (*this); }
	
	int size() const { return directProduct(bounds); }
	
	/** If some bound less or equal to zero, all next operations are invalid */
	void assertBoundsValid() const {
		for (int i = 0; i < D; i++) {
			assert_gt(bounds(i), 0);
		}
	}
	
protected:
	const IntD bounds;
	
	/**
	 * Increase index at position i.
	 * @return is bound at position i reached
	 */
	bool increment(const int i) {
		(*this)(i) = ((*this)(i) + 1) % this->bounds(i);
		return !(bool)((*this)(i));
	}
};


/**
 * Multiindex where first index is "slow" and last is "fast",
 * i.e. ++SlowXFastZ({0, 0}, {5, 5}) == SlowXFastZ({0, 1}, {5, 5}).
 * Performs iteration from Zeros() to sizes.
 */
template<int D> struct SlowXFastZ;

template<>
struct SlowXFastZ<1> : public Multiindex<1> {
	typedef Multiindex<1> Base;
	typedef Base::IntD    IntD;
	using Base::Base;
	
	SlowXFastZ& operator++() {
		++((*this)(0));
		return (*this);
	}

	static SlowXFastZ begin(const IntD& bounds_) {
		return SlowXFastZ({0}, bounds_);
	}
	static SlowXFastZ end(const IntD& bounds_) {
		return SlowXFastZ(bounds_, bounds_);
	}
	
	SlowXFastZ begin() const { return begin(this->bounds); }
	SlowXFastZ end() const { return end(this->bounds); }
};

template<>
struct SlowXFastZ<2> : public Multiindex<2> {
	typedef Multiindex<2> Base;
	typedef Base::IntD    IntD;
	using Base::Base;
	
	SlowXFastZ& operator++() {
		this->increment(1) && (++((*this)(0)));
		return (*this);
	}

	static SlowXFastZ begin(const IntD& bounds_) {
		return SlowXFastZ({0, 0}, bounds_);
	}
	static SlowXFastZ end(const IntD& bounds_) {
		return SlowXFastZ({bounds_(0), 0}, bounds_);
	}
	
	SlowXFastZ begin() const { return begin(this->bounds); }
	SlowXFastZ end() const { return end(this->bounds); }
};

template<>
struct SlowXFastZ<3> : public Multiindex<3> {
	typedef Multiindex<3> Base;
	typedef Base::IntD    IntD;
	using Base::Base;
	
	SlowXFastZ& operator++() {
		this->increment(2) && this->increment(1) && ++((*this)(0));
		return (*this);
	}

	static SlowXFastZ begin(const IntD& bounds_) {
		return SlowXFastZ({0, 0, 0}, bounds_);
	}
	static SlowXFastZ end(const IntD& bounds_) {
		return SlowXFastZ({bounds_(0), 0, 0}, bounds_);
	}
	
	SlowXFastZ begin() const { return begin(this->bounds); }
	SlowXFastZ end() const { return end(this->bounds); }
};


/**
 * Multiindex where first index is "fast" and last is "slow",
 * i.e. ++SlowZFastX({0, 0}, {5, 5}) == SlowXFastZ({1, 0}, {5, 5}).
 * Performs iteration from Zeros() to sizes.
 */
template<int D> struct SlowZFastX;

template<>
struct SlowZFastX<1> : public Multiindex<1> {
	typedef Multiindex<1> Base;
	typedef Base::IntD    IntD;
	using Base::Base;
	
	SlowZFastX& operator++() {
		++((*this)(0));
		return (*this);
	}

	static SlowZFastX begin(const IntD& bounds_) {
		return SlowZFastX({0}, bounds_);
	}
	static SlowZFastX end(const IntD& bounds_) {
		return SlowZFastX(bounds_, bounds_);
	}
	
	SlowZFastX begin() const { return begin(this->bounds); }
	SlowZFastX end() const { return end(this->bounds); }
};


template<>
struct SlowZFastX<2> : public Multiindex<2> {
	typedef Multiindex<2> Base;
	typedef Base::IntD    IntD;
	using Base::Base;
	
	SlowZFastX& operator++() {
		this->increment(0) && (++((*this)(1)));
		return (*this);
	}

	static SlowZFastX begin(const IntD& bounds_) {
		return SlowZFastX({0, 0}, bounds_);
	}
	static SlowZFastX end(const IntD& bounds_) {
		return SlowZFastX({0, bounds_(1)}, bounds_);
	}
	
	SlowZFastX begin() const { return begin(this->bounds); }
	SlowZFastX end() const { return end(this->bounds); }
};


template<>
struct SlowZFastX<3> : public Multiindex<3> {
	typedef Multiindex<3> Base;
	typedef Base::IntD    IntD;
	using Base::Base;
	
	SlowZFastX& operator++() {
		this->increment(0) && this->increment(1) && (++((*this)(2)));
		return (*this);
	}

	static SlowZFastX begin(const IntD& bounds_) {
		return SlowZFastX({0, 0, 0}, bounds_);
	}
	static SlowZFastX end(const IntD& bounds_) {
		return SlowZFastX({0, 0, bounds_(2)}, bounds_);
	}
	
	SlowZFastX begin() const { return begin(this->bounds); }
	SlowZFastX end() const { return end(this->bounds); }
};


/**
 * Iteration through a box (cube) from min to max
 * @note if min == Zeros(), iterator is same to its RelativeIteratorType<D>
 */
template<int D, template<int> class RelativeIteratorType>
struct BoxIterator : public VectorInt<D> {
	typedef VectorInt<D>               IntD;
	typedef RelativeIteratorType<D>    RelativeIterator;
	
	/**
	 * Iterator from min INclusive to max EXclusive.
	 * Initial position at start, start MUST be between min and max.
	 */
	BoxIterator(const IntD& start, const IntD& min, const IntD& max) :
			IntD(start), shift(min), relativeIterator(start - min, max - min) { }
	
	BoxIterator(const RelativeIteratorType<D>& relIter, const IntD shift_) :
			IntD(relIter + shift_), shift(shift_), relativeIterator(relIter) { }
	
	
	BoxIterator& operator++() {
		++relativeIterator;
		this->copyFrom(relativeIterator + shift);
		return (*this);
	}
	
	
	BoxIterator begin() const {
		return BoxIterator(relativeIterator.begin(), shift);
	}
	BoxIterator end() const {
		return BoxIterator(relativeIterator.end(), shift);
	}
	
	
	int size() const { return relativeIterator.size(); }
	
	void assertBoundsValid() const {
		relativeIterator.assertBoundsValid();
	}
	
protected:
	const IntD shift;
	RelativeIterator relativeIterator;
	
};


///@}

}
}

#endif // LIBGCM_LINAL_MULTIINDEX_HPP
