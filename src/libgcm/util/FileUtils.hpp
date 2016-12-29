#ifndef LIBGCM_FILEUTILS_HPP
#define LIBGCM_FILEUTILS_HPP


#include <fstream>
#include <vector>

#include <libgcm/util/infrastructure/infrastructure.hpp>

namespace gcm {
class FileUtils {
public:
// FIXME clean this file as well as snapshotters
	
	template<typename T>
	static void writeToTextFile(const std::string& fileName, const T& t) {
		std::ofstream fileStream(fileName, std::ios::out);
		assert_true(fileStream.is_open());
		fileStream << t;
		assert_true(fileStream.good());
		fileStream.close();
	}
	
	template<typename T>
	static void writeStdVectorsToTextFile(const std::string fileName,
			const std::vector<std::vector<T>>& vecs) {
		std::ofstream fileStream;
		openTextFileStream(fileStream, fileName);
		writeStdVectorsToTextFileStream(fileStream, vecs);
		closeFileStream(fileStream);
	}
	

	static void openBinaryFileStream(std::ofstream& fileStream,
	                                 const std::string& fileName) {
		fileStream.open(fileName, std::ios::binary);
		assert_true(fileStream.is_open());
		fileStream.clear();
	}

	
	/** Open file stream to write in */
	static void openTextFileStream(std::ofstream& fileStream,
	                               const std::string& fileName) {
		fileStream.open(fileName, std::ios::out);
		assert_true(fileStream.is_open());
		fileStream.clear();
	}
	
	
	/** Open file stream to read from */
	static void openTextFileStream(std::ifstream& fileStream,
	                               const std::string& fileName) {
		fileStream.open(fileName, std::ios::in);
		assert_true(fileStream.is_open());
	}


	template<typename T>
	static void writeStdVectorToTextFileStream(std::ofstream& fileStream,
	                                           const std::vector<T>& vec) {
		for (const auto& it : vec) {
			fileStream << it << std::endl;
		}
	}
	
	
	template<typename T>
	static void writeStdVectorsToTextFileStream(std::ofstream& fileStream,
			const std::vector<std::vector<T>>& vecs) {
		assert_false(vecs.empty());
		size_t s = vecs[0].size();
		for (const auto& v : vecs) {
			assert_eq(s, v.size());
		}
		for (size_t i = 0; i < s; i++) {
			for (const auto& v : vecs) {
				fileStream << v[i] << "\t";
			}
			fileStream << std::endl;
		}
	}


	template<typename T>
	static void writeStdVectorToBinaryFileStream(std::ofstream& fileStream,
	                                             const std::vector<T>& vec) {
		writeArrayToBinaryFileStream(fileStream, &(vec[0]), vec.size());
	}


	template<typename T>
	static void writeArrayToBinaryFileStream(std::ofstream& fileStream,
	                                         const T* array, const size_t sizeOfArray) {
		assert_gt(sizeOfArray, 0);
		auto bufferSize = (std::streamsize)(sizeOfArray * sizeof(T));
		auto previousNumberOfBytes = fileStream.tellp();
		assert_true(fileStream.write(reinterpret_cast<const char*>(array), bufferSize));
		auto currentNumberOfBytes = fileStream.tellp();
		assert_eq(bufferSize, currentNumberOfBytes - previousNumberOfBytes);
	}


	template<typename PlaceToWriteInput>
	static void readFromTextFile(const std::string& fileName,
			PlaceToWriteInput& placeToWriteInput) {
		
		std::ifstream inputFileStream(fileName);
		assert_true(inputFileStream.is_open());
		
		inputFileStream >> placeToWriteInput;
		closeFileStream(inputFileStream);
	}


	static void closeFileStream(std::ifstream& fileStream) {
		assert_true(fileStream.good());
		fileStream.close();
	}
	
	
	static void closeFileStream(std::ofstream& fileStream) {
		assert_true(fileStream.good());
		fileStream.close();
	}
	
};


}

#endif /* LIBGCM_FILEUTILS_HPP */
