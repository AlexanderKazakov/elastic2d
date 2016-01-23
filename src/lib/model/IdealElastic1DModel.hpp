#ifndef LIBGCM_IDEALELASTIC1DMODEL_HPP
#define LIBGCM_IDEALELASTIC1DMODEL_HPP

#include <lib/model/Model.hpp>
#include <lib/nodes/Node.hpp>

namespace gcm {
	class IdealElastic1DModel : public Model {
	public:
		typedef IdealElastic1DNode Node;
		typedef GcmMatrices<2, 1, IsotropicMaterial> GCM_MATRICES;
		static const int DIMENSIONALITY = GCM_MATRICES::DIMENSIONALITY;

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

		static real GetVx(const Node& node) { return node.u.V[0]; };
		static real GetSxx(const Node& node) { return node.u.S[0]; };
		static real GetPressure(const Node& node) { return node.u.getPressure(); };

		static void SetVx(const real& value, Node& node) { node.u.V[0] = value; };
		static void SetSxx(const real& value, Node& node) { node.u.S[0] = value; };
		static void SetPressure(const real& value, Node& node) { node.u.setPressure(value); };

	};
}


#endif // LIBGCM_IDEALELASTIC1DMODEL_HPP
