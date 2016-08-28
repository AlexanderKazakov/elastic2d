#ifndef LIBGCM_CUBIC_BORDERCONDITIONS_HPP
#define LIBGCM_CUBIC_BORDERCONDITIONS_HPP

#include <libgcm/util/math/Area.hpp>

namespace gcm {
namespace cubic {


/**
 * Applying border conditions for cubic meshes.
 * The approach is ghost (fixture) nodes on borders.
 * Setting appropriate values in ghost nodes before calculation,
 * one can calculate borders in the same manner as inner nodes.
 */
template<typename Mesh>
struct BorderConditions {
	typedef typename Mesh::PdeVector                        PdeVector;
	typedef typename Mesh::Iterator                         Iterator;
	typedef typename Mesh::PdeVariables                     PdeVariables;
	
	typedef std::function<real(real)>                       TimeDependency;
	typedef std::map<PhysicalQuantities::T, TimeDependency> Map;
	
	struct Condition {
		Condition(const std::shared_ptr<Area> area_, const Map& values_) :
				area(area_), values(values_) {
			for (const auto& q : values) {
				assert_eq(PdeVariables::QUANTITIES.count(q.first), 1);
			}
		}
		std::shared_ptr<Area> area; ///< area of surface
		Map values; ///< fixed values
	};
	
	
	BorderConditions(const Task& task) {
		for (const auto& bc : task.cubicGridBorderConditions) {
			conditions.push_back(Condition(bc.area, bc.values));
		}
	}
	
	
	void apply(Mesh* mesh, const int stage) const {
		/// handling borders
		if (stage == 0) {
		// special for x-axis (because MPI partitioning along x-axis)
			if (Mpi::ForceSequence() || Mpi::Rank() == 0) {
				handleSide(mesh, stage, false);
			}
			if (Mpi::ForceSequence() || Mpi::Rank() == Mpi::Size() - 1) {
				handleSide(mesh, stage, true);
			}
			return;
		}
		// for other axes
		handleSide(mesh, stage, false);
		handleSide(mesh, stage, true);
	}
	
	void handleSide(
			Mesh* mesh, const int direction, const bool onTheRight) const {
		
		auto borderIter = onTheRight ?
				mesh->slice(direction, mesh->sizes(direction) - 1) :
				mesh->slice(direction, 0);
		
		while (borderIter != borderIter.end()) {
			for (const auto& condition : conditions) {
				if (condition.area->contains(mesh->coords(borderIter))) {
					handleBorderPoint(mesh, borderIter, condition.values, direction, onTheRight);
				}
			}
			++borderIter;
		}
	}
	
	static void handleBorderPoint(Mesh* mesh, const Iterator& borderIter,
			const Map& values, const int direction, const bool onTheRight) {
		
		int innerSign = onTheRight ? -1 : 1;
		for (int a = 1; a <= mesh->borderSize; a++) {
			auto realIter = borderIter; realIter(direction) += innerSign * a;
			auto virtIter = borderIter; virtIter(direction) -= innerSign * a;
	
			mesh->_pde(virtIter) = mesh->pde(realIter);
			for (const auto& q : values) {
				const auto& quantity = q.first;
				const auto& timeDependency = q.second;
				
				real realValue = PdeVariables::QUANTITIES.at(quantity).
						Get(mesh->pde(realIter));
				real virtValue = -realValue + 2 * timeDependency(Clock::Time());
				
				PdeVariables::QUANTITIES.at(quantity).
						Set(virtValue, mesh->_pde(virtIter));
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
