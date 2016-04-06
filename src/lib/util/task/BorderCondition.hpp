#ifndef LIBGCM_BORDERCONDITIONS_SIGMA_HPP
#define LIBGCM_BORDERCONDITIONS_SIGMA_HPP

#include <lib/linal/linal.hpp>
#include <lib/mesh/DefaultMesh.hpp>
#include <lib/Engine.hpp>


namespace gcm {

/**
 * Border condition.
 * Ideas and designations from Chelnokov's PhD thesis.
 * 
 * General expression of linear border conditions is
 * \f$ B \vec{u}_{n+1} = \vec{b}(t_{n+1}) \f$.
 * 
 * Value of vector \vec{b} is kept in local basis \f$ {\vec{\tau}, \vec{n}} \f$,
 * connected with body border. More about the basis @see linal::createLocalBasis
 * 
 * I.e, to set normal velocity out of the body,
 * let \f$ \vec{b} = (0, velocity) \f$.
 * To set normal force into the body,
 * let \f$ \vec{b} = (0, -force) \f$.
 * 
 * It can be done more generally: not for tasks like
 * "fixed_force" but for any list of quantities with their time dependencies,
 * but it requires some logic about basis in matrix B (TODO)
 */
template<typename TModel>
class BorderCondition {
public:
	static const int PDE_SIZE = TModel::PDE_SIZE;
	static const int DIMENSIONALITY = TModel::DIMENSIONALITY;
	/// Legal number of outer characteristics in border node
	static const int OUTER_NUMBER = DIMENSIONALITY;
	
	typedef typename TModel::PdeVector            PdeVector;
	typedef linal::Matrix<OUTER_NUMBER, PDE_SIZE> BorderMatrix;
	typedef linal::Vector<OUTER_NUMBER>           BorderVector;
	typedef std::function<BorderVector(real)>     BorderVectorTimeDependency;
	
	BorderCondition(const Statement::BorderCondition& task) :
			type_(task.type), b_(createBorderVectorTimeDependency(task.values)) {
		switch (type_) {
			case gcm::BorderConditions::T::FIXED_VELOCITY:
			case gcm::BorderConditions::T::FIXED_FORCE:
				break;
			default:
				THROW_INVALID_ARG("Unknown type of border condition");
		}
	}
	
	/**
	 * General expression of linear border conditions is
	 * \f$ B \vec{u}_{n+1} = \vec{b}(t_{n+1}) \f$
	 */
	BorderMatrix B(const linal::Vector<DIMENSIONALITY>& normal) const {
		switch (type_) {
			case gcm::BorderConditions::T::FIXED_VELOCITY:
				return borderMatrixFixedVelocity(normal);
				break;
			case gcm::BorderConditions::T::FIXED_FORCE:
				return borderMatrixFixedForce(normal);
				break;
			default:
				THROW_INVALID_ARG("This is impossible");
		}
	}
	
	/**
	 * General expression of linear border conditions is
	 * \f$ B \vec{u}_{n+1} = \vec{b}(t_{n+1}) \f$
	 */
	BorderVector b(const linal::Vector<DIMENSIONALITY>& normal) const {
		auto S = linal::createLocalBasis(normal); // transfer matrix
		return S * b_(Clock::TimeAtNextTimeLayer());
	}
	
private:
	const gcm::BorderConditions::T type_;
	const BorderVectorTimeDependency b_;
	
	BorderMatrix borderMatrixFixedForce
	(const linal::Vector<DIMENSIONALITY>& normal) const {
	///	Set border matrix for the case of fixed force on border
		BorderMatrix B_;
		for (int i = 0; i < DIMENSIONALITY; i++) {
			PdeVector pde = PdeVector::zeros();
			for (int j = 0; j < DIMENSIONALITY; j++) {
				pde.sigma(i, j) = normal(j);
			}
			B_.setString(i, pde);
		}
		return B_;
	}
	
	BorderMatrix borderMatrixFixedVelocity
	(const linal::Vector<DIMENSIONALITY>&) const {
	///	Set border matrix for the case of fixed velocity on border.
	/// Velocity border matrix B is independent from border normal
		BorderMatrix B_;
		for (int i = 0; i < DIMENSIONALITY; i++) {
			PdeVector pde = PdeVector::zeros();
			pde.V[i] = 1;
			B_.setString(i, pde);
		}
		return B_;
	}
	
	BorderVectorTimeDependency createBorderVectorTimeDependency
	(const std::vector<Statement::TimeDependency> values) const {
	/// translate border condition from task format to own format
		
		assert_eq(values.size(), OUTER_NUMBER);
		Statement::TimeDependency v[OUTER_NUMBER];
		for (int i = 0; i < OUTER_NUMBER; i++) {
			v[i] = values[(size_t)i];
		}
		
		BorderVectorTimeDependency ans = [v](real t) {
			BorderVector borderVector;
			for (int i = 0; i < OUTER_NUMBER; i++) {
				borderVector(i) = v[i](t);
			}
			return borderVector;
		};
		
		return ans;
	}
	
};


}

#endif // LIBGCM_BORDERCONDITIONS_SIGMA_HPP
