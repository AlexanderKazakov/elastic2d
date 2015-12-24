#ifndef LIBGCM_NODE_HPP
#define LIBGCM_NODE_HPP

#include <memory>

#include "lib/PDEMatrices.hpp"

namespace gcm {
	class Node {
	public:
		Vector u;
		std::shared_ptr <PDEMatrices> matrix;

		inline const real &get(const NodeMap nodeMap) const {
			return u.get(static_cast<int>(nodeMap));
		};

		inline real &operator()(const NodeMap nodeMap) {
			return u(static_cast<int>(nodeMap));
		};
	};
}

#endif //LIBGCM_NODE_HPP
