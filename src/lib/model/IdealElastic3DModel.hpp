#ifndef LIBGCM_IDEALELASTIC3DMODEL_HPP
#define LIBGCM_IDEALELASTIC3DMODEL_HPP

#include "lib/model/Model.hpp"
#include "lib/nodes/IdealElastic3DNode.hpp"
#include "lib/gcm_matrices/IsotropicElastic3DGcmMatrices.hpp"

namespace gcm {
	class IdealElastic3DModel : public Model {
	public:
		typedef IdealElastic3DNode Node;
		typedef IsotropicElastic3DGcmMatrices GcmMatrices;
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
		static real GetVy(const Node& node) { return node.u.Vy; };
		static real GetVz(const Node& node) { return node.u.Vz; };
		static real GetSxx(const Node& node) { return node.u.Sxx; };
		static real GetSxy(const Node& node) { return node.u.Sxy; };
		static real GetSxz(const Node& node) { return node.u.Sxz; };
		static real GetSyy(const Node& node) { return node.u.Syy; };
		static real GetSyz(const Node& node) { return node.u.Syz; };
		static real GetSzz(const Node& node) { return node.u.Szz; };
		static real GetPressure(const Node& node) { return node.getPressure(); };

		static void SetVx(const real& value, Node& node) { node.u.Vx = value; };
		static void SetVy(const real& value, Node& node) { node.u.Vy = value; };
		static void SetVz(const real& value, Node& node) { node.u.Vz = value; };
		static void SetSxx(const real& value, Node& node) { node.u.Sxx = value; };
		static void SetSxy(const real& value, Node& node) { node.u.Sxy = value; };
		static void SetSxz(const real& value, Node& node) { node.u.Sxz = value; };
		static void SetSyy(const real& value, Node& node) { node.u.Syy = value; };
		static void SetSyz(const real& value, Node& node) { node.u.Syz = value; };
		static void SetSzz(const real& value, Node& node) { node.u.Szz = value; };
		static void SetPressure(const real& value, Node& node) { node.setPressure(value); };

	};
}


#endif // LIBGCM_IDEALELASTIC3DMODEL_HPP
