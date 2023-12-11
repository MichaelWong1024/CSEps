# Compiler settings
CC = gcc
CFLAGS = -Wall -Wextra -std=c99

# Libraries to link
LIBS = -lm

# Source file
SOURCE = newfloat.c

# Executable name
EXECUTABLE = newfloat

# Default rule
all: $(EXECUTABLE)

# Linking the program
$(EXECUTABLE): $(SOURCE)
	$(CC) $(CFLAGS) $(SOURCE) -o $(EXECUTABLE) $(LIBS)

# Clean up
clean:
	rm -f $(EXECUTABLE)

