# Define the compiler
CC=gcc

# Define any compile-time flags
CFLAGS=-Wall -fopenmp

# Define the source file and the executable name
SOURCE=bellmanford.c
EXECUTABLE=bellmanford

# Default target
all: $(EXECUTABLE)

# Rule for building the executable
$(EXECUTABLE): $(SOURCE)
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(SOURCE)

# Rule for cleaning up
clean:
	rm -f $(EXECUTABLE)

