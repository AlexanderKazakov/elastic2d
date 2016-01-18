#ifndef LIBGCM_IDEALELASTIC1DMODEL_HPP
#define LIBGCM_IDEALELASTIC1DMODEL_HPP

#include "lib/model/Model.hpp"
#include "lib/nodes/IdealElastic1DNode.hpp"
#include "lib/gcm_matrices/Elastic1DGcmMatrices.hpp"

namespace gcm {
	class IdealElastic1DModel : public Model {
	public:
		typedef IdealElastic1DNode Node;
		typedef Elastic1DGcmMatrices GcmMatrices;
		static const int DIMENSIONALITY = GcmMatrices::DIMENSIONALITY;

		// The code below is to provide for classes like snapshotters, initial condition setters, etc
		// unified for all models interface to information about model's physical quantities
		// and access to that values in some node.
		// The function pointers are supposed to reduce call overhead.
		// The map shouldn't be used at every access to every node, but just once before handling
		// a large portion of nodes
		typedef real (*Getter)(const Node& node);
		typedef void (*Setter)(const real& value, Node& node);
		struct GetSet {
			GetSet(Getter Get, Setter Set) : Get(Get), Set(Set) { };
			Getter Get;
			Setter Set;
		};
		static const std::map<PhysicalQuantities::T, GetSet> QUANTITIES;

		static real GetVx(const Node& node) { return node.u.Vx; };
		static real GetSxx(const Node& node) { return node.u.Sxx; };
		static real GetPressure(const Node& node) { return node.getPressure(); };

		static void SetVx(const real& value, Node& node) { node.u.Vx = value; };
		static void SetSxx(const real& value, Node& node) { node.u.Sxx = value; };
		static void SetPressure(const real& value, Node& node) { node.setPressure(value); };

	};
}


#endif // LIBGCM_IDEALELASTIC1DMODEL_HPP
