COMPILER=mpic++
CFLAGS=-O3 --std=c++11 -Isrc
LIBRARIES=-L/usr/lib
LIBFLAGS=-lm

EXECUTABLE=build/elastic2d
SOURCES=$(shell find src/lib -name '*.cpp') src/main.cpp
OBJECTS=$(SOURCES:%.cpp=%.o)

all: $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS)
	$(COMPILER) $(CFLAGS) $(INCLUDES) $(LIBRARIES) $(OBJECTS) -o $@ $(LIBFLAGS)
	
%.o: %.cpp
	$(COMPILER) $(CFLAGS) $(INCLUDES) $(LIBRARIES) -c $< -o $@ $(LIBFLAGS)

clean:
	rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)
