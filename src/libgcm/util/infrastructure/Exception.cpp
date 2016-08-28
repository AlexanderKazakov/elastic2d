#include <libgcm/util/infrastructure/Exception.hpp>

using namespace gcm;

using std::string;

const string& Exception::
getCallStack() const {
	return callStack;
}


Exception::
Exception(int _code, const string& _message, const string& _file, int _line) {
	this->code = _code;
	this->line = _line;
	this->message = _message;
	this->file = _file;

	// save call stack
	void** buffer = new void*[100];
	int callStackLen = backtrace(buffer, 100);
	char** cs = backtrace_symbols(buffer, callStackLen);


	auto demangle = [] (char* symbol)->string
	{
		size_t size;
		int status;
		char temp[128];
		char* demangled;

		string result;

		// first, try to demangle a c++ name
		if (1 == sscanf(symbol, "%*[^(]%*[^_]%127[^)+]", temp)) {
			if ((demangled = abi::__cxa_demangle(temp, NULL, &size, &status))) {
				result = string(demangled);
				free(demangled);
			}
		} else if (1 == sscanf(symbol, "%127s", temp)) {
			// if that didn't work, try to get a regular c symbol
			result = string(temp);
		} else {
			// if all else fails, just return the symbol
			result = string(symbol);
		}
		return result;
	};

	callStack = "";
	for (int i = 1; i < callStackLen; i++) {
		callStack += demangle(cs[i]) + "\n";
	}

	free(cs);
	delete[] buffer;
}


int Exception::
getCode() const {
	return code;
}


const string& Exception::
getMessage() const {
	return message;
}


const string& Exception::
getFile() const {
	return file;
}


int Exception::
getLine() const {
	return line;
}


const std::string Exception::
what() const {
	return "Exception was thrown: " + getMessage() + "\n @" + getFile() + ":" + std::to_string(
	               getLine()) + "\nCall stack: \n" + getCallStack();
}


