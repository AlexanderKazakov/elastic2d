#ifndef LIBGCM_BORDERCONDITION_HPP
#define LIBGCM_BORDERCONDITION_HPP

#include <libgcm/linal/linal.hpp>
#include <libgcm/util/task/Task.hpp>
#include <libgcm/engine/GlobalVariables.hpp>


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
	typedef typename TModel::BorderMatrix         BorderMatrix;
	typedef typename TModel::BorderVector         BorderVector;
	typedef std::function<BorderVector(real)>     BorderVectorTimeDependency;
	
	const BorderConditions::T type;
	
	BorderCondition(const Task::BorderCondition& task) :
			type(task.type),
			b_(createBorderVectorTimeDependency(task.values)) { }
	
	/**
	 * General expression of linear border conditions is
	 * \f$ B \vec{u}_{n+1} = \vec{b}(t_{n+1}) \f$.
	 * @param time time at (n+1)'th layer (not at n'th!)
	 * @return \vec{b} in local (connected with border) basis
	 */
	BorderVector b(const real time) const {
		return b_(time);
	}
	
	
private:
	const BorderVectorTimeDependency b_;
	
	
	BorderVectorTimeDependency createBorderVectorTimeDependency
	(const std::vector<Task::TimeDependency> values) const {
	/// translate border condition from task format to own format
		
		assert_eq(values.size(), OUTER_NUMBER);
		Task::TimeDependency v[OUTER_NUMBER];
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

#endif // LIBGCM_BORDERCONDITION_HPP
