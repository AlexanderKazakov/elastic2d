#ifndef LIBGCM_CUBIC_BORDERCONDITIONS_HPP
#define LIBGCM_CUBIC_BORDERCONDITIONS_HPP

#include <libgcm/util/task/Task.hpp>
#include <libgcm/grid/AbstractGrid.hpp>


namespace gcm {
namespace cubic {

/**
 * Applying border conditions for cubic meshes.
 * The approach is ghost (fixture) nodes on borders.
 * Setting appropriate values in ghost nodes before calculation,
 * we can calculate borders in the same manner as inner nodes.
 */
class AbstractBorderConditions {
public:
	virtual void apply(AbstractGrid& mesh_, const int direction) const = 0;
};


template<typename Mesh>
class BorderConditions : public AbstractBorderConditions {
public:
	typedef typename Mesh::PdeVector                        PdeVector;
	typedef typename Mesh::Iterator                         Iterator;
	typedef typename Mesh::PdeVariables                     PdeVariables;
	typedef typename Mesh::Grid::PartIterator               PartIterator;
	
	typedef std::function<real(real)>                       TimeDependency;
	typedef std::map<PhysicalQuantities::T, TimeDependency> Map;
	
	struct Condition {
		/// direction of border normal (x=0, y=1, z=2) (aka stage);
		/// before a concrete stage, we have to apply border condition
		/// along this direction only (others are unnecessary)
		int direction;
		/// lists of nodes to perform border condition on
		std::vector<Iterator> leftNodes, rightNodes;
		/// list of physical quantities in border condition
		Map values;
	};
	
	
	BorderConditions(const Task& task, const AbstractGrid& mesh_) {
		const Mesh& mesh = dynamic_cast<const Mesh&>(mesh_);
		const auto meshConditions = task.cubicBorderConditions.find(mesh.id);
		if (meshConditions == task.cubicBorderConditions.end()) { return; }
		
		for (const Task::CubicBorderCondition& bc : meshConditions->second) {
			Condition condition;
			condition.direction = bc.direction;
			condition.values = bc.values;
			
			/// check specified physical quantities presence in pde variables
			for (const auto quantity : bc.values) {
				assert_false(PdeVariables::QUANTITIES.find(quantity.first)
						== PdeVariables::QUANTITIES.end());
			}
			
			/// find border nodes to apply conditions
			for (PartIterator leftNode = mesh.leftBorder(condition.direction);
					leftNode  != leftNode .end(); ++leftNode) {
				if (bc.area->contains(mesh.coords(leftNode))) {
					condition.leftNodes.push_back(leftNode);
				}
			}
			for (PartIterator rightNode = mesh.rightBorder(condition.direction);
					rightNode != rightNode.end(); ++rightNode) {
				if (bc.area->contains(mesh.coords(rightNode))) {
					condition.rightNodes.push_back(rightNode);
				}
			}
			
			conditions.push_back(condition);
		}
	}
	
	
	virtual void apply(AbstractGrid& mesh_, const int direction) const override {
		Mesh& mesh = dynamic_cast<Mesh&>(mesh_);
		for (const Condition& c : conditions) if (c.direction == direction) {
			for (const Iterator& left : c.leftNodes) {
				handleBorderPoint(mesh, left, c.values, direction, 1);
			}
			for (const Iterator& right : c.rightNodes) {
				handleBorderPoint(mesh, right, c.values, direction, -1);
			}
		}
	}
	
	
	static void handleBorderPoint(Mesh& mesh, const Iterator& iter,
			const Map& values, const int direction, const int innerSign) {
		
		for (int a = 1; a <= mesh.borderSize; a++) {
			auto inner = iter; inner(direction) += innerSign * a;
			auto ghost = iter; ghost(direction) -= innerSign * a;
			
			mesh._pde(ghost) = mesh.pde(inner);
			for (const auto& q : values) {
				const auto& quantity = q.first;
				const auto& timeDependency = q.second;
				
				real innerValue = PdeVariables::QUANTITIES.at(quantity).
						Get(mesh.pde(inner));
				real ghostValue = -innerValue + 2 * timeDependency(Clock::Time());
				
				PdeVariables::QUANTITIES.at(quantity).
						Set(ghostValue, mesh._pde(ghost));
			}
		}
	}
	
	
private:
	/// list of border conditions applied in sequence (overwriting previous)
	std::vector<Condition> conditions;
	
};


} // namespace cubic
} // namespace gcm


#endif // LIBGCM_CUBIC_BORDERCONDITIONS_HPP
