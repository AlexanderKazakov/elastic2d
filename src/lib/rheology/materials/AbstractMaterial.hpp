#ifndef LIBGCM_ABSTRACTMATERIAL_HPP
#define LIBGCM_ABSTRACTMATERIAL_HPP

#include <initializer_list>

#include <lib/util/Enum.hpp>
#include <lib/linal/linal.hpp>

namespace gcm {
	struct AbstractMaterial {
		virtual ~AbstractMaterial() = default;
	};
}

#endif // LIBGCM_ABSTRACTMATERIAL_HPP
