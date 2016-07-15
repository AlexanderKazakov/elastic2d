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
 */
template<typename TModel>
class BorderCondition {
public:
	static const int PDE_SIZE = TModel::PDE_SIZE;
	static const int DIMENSIONALITY = TModel::DIMENSIONALITY;
	/// Number of characteristics with slopes of the same sign.
	/// It is equal to number of outer characteristics in border node.
	static const int OUTER_NUMBER = TModel::OUTER_NUMBER;
	
	typedef typename TModel::PdeVariables         PdeVariables;
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
	BorderMatrix B(const linal::Vector<DIMENSIONALITY>& borderNormal) const {
		switch (type_) {
			case gcm::BorderConditions::T::FIXED_VELOCITY:
				return TModel::borderMatrixFixedVelocity(borderNormal);
				break;
			case gcm::BorderConditions::T::FIXED_FORCE:
				return TModel::borderMatrixFixedForce(borderNormal);
				break;
			default:
				THROW_INVALID_ARG("Unknown type of border condition");
		}
	}
	
	
	/**
	 * General expression of linear border conditions is
	 * \f$ B \vec{u}_{n+1} = \vec{b}(t_{n+1}) \f$.
	 * @return \vec{b} in local (connected with border) basis
	 */
	BorderVector b() const {
		return b_(Clock::TimeAtNextTimeLayer());
	}
	
	
private:
	const gcm::BorderConditions::T type_;
	const BorderVectorTimeDependency b_;
	
	
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
