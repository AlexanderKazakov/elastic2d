#ifndef LIBGCM_STRINGUTILS_HPP
#define LIBGCM_STRINGUTILS_HPP

#include <iomanip>
#include <sstream>

namespace gcm {
	class StringUtils {
	public:
		/**
		 * Print number to string of specified length with leading zeroes
		 */
		static std::string toString(const int number, const int length) {
			std::ostringstream strStream;
			strStream << std::setfill('0') << std::setw(length) << number;
			return strStream.str();
		};
	};
}

#endif /* LIBGCM_STRINGUTILS_HPP */
