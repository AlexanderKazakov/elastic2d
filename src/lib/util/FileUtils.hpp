#ifndef LIBGCM_FILEUTILS_HPP
#define LIBGCM_FILEUTILS_HPP

#include <fstream>
#include <vector>

#include <lib/util/Assertion.hpp>

namespace gcm {
	class FileUtils {
	public:
		static void openBinaryFileStream(std::ofstream& fileStream, 
		                                 const std::string& fileName) {
			fileStream.open(fileName, std::ios::binary);
			assert_true(fileStream.is_open());
			fileStream.clear();
		};
		static void openTextFileStream(std::ofstream& fileStream, 
		                               const std::string& fileName) {
			fileStream.open(fileName, std::ios::out);
			assert_true(fileStream.is_open());
			fileStream.clear();
		};
		
		template<typename T>
		static void writeStdVectorToTextFileStream(std::ofstream& fileStream,
		                                           const std::vector<T>& vec) {
			for (const auto& t : vec) {
				fileStream << t << std::endl;
			}
		};
		
		template<typename T>
		static void writeStdVectorToBinaryFileStream(std::ofstream& fileStream,
		                                             const std::vector<T>& vec) {
			writeArrayToBinaryFileStream(fileStream, &(vec[0]), vec.size());
		};
		template<typename T>
		static void writeArrayToBinaryFileStream(std::ofstream& fileStream,
				const T* array, const size_t sizeOfArray) {
			assert_gt(sizeOfArray, 0);
			auto bufferSize = (std::streamsize)(sizeOfArray * sizeof(T));
			auto previousNumberOfBytes = fileStream.tellp();
			assert_true(fileStream.write(reinterpret_cast<const char*>(array), bufferSize));
			auto currentNumberOfBytes = fileStream.tellp();
			assert_eq(bufferSize, currentNumberOfBytes - previousNumberOfBytes);
		};
		
		static void closeFileStream(std::ofstream& fileStream) {
			assert_true(fileStream.good());
			fileStream.close();
		};
	};
}

#endif /* LIBGCM_FILEUTILS_HPP */
