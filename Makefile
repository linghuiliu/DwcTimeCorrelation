TARGET = DwcTimeCorrelation
LIBS = -lm
CC = g++
CFLAGS = -Wall -O3
ROOTFLAG = `root-config --cflags`
ROOTLIBS = `root-config --libs`
#INCLUDE = -I /cern/root-v5.34.28-gcc44/include

.PHONY: clean all default

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
HEADERS = $(wildcard *.h)

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) $(ROOTFLAG) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(ROOTLIBS) $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f $(TARGET)
