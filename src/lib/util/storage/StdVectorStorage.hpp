#ifndef LIBGCM_STDVECTORSTORAGE_HPP
#define LIBGCM_STDVECTORSTORAGE_HPP

#include <vector>

namespace gcm {
	template<typename T>
	struct StdVectorStorage : public std::vector<T> {
		void zeroInitialize(size_t InitSize) {
			this->resize(InitSize);
			memset(&((*this)[0]), 0, this->size() * sizeof(T));
		};
	};
}

#endif // LIBGCM_STDVECTORSTORAGE_HPP
