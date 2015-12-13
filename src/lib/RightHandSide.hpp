#ifndef ELASTIC2D_RIGHTHANDSIDE_HPP
#define ELASTIC2D_RIGHTHANDSIDE_HPP

#include "lib/Node.hpp"

/**
 * The PDE is:
 * d\vec{u}/dt + \mat{Ax} \partial{\vec{u}}/\partial{x} + \mat{Ay} \partial{\vec{u}}/\partial{y} = \vec{f}.
 * The class represents the \vec{f}
 */
class RightHandSide {
public:
	/**
	 * @return the value of right-hand side of PDE in @node at time @t
	 */
	Vector getRightHandSide(const Node& node, const double& t) const;
};


#endif //ELASTIC2D_RIGHTHANDSIDE_HPP
