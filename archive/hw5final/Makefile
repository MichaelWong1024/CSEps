# Makefile for compiling bellmanford program

# Compiler to use
CC = gcc

# Compiler flags
CFLAGS = -Wall -g -fopenmp

# Target executable name
TARGET = bellmanford

# Object files
OBJS = main.o bellmanford.o

# Rule for making the target
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

# Rule for making main.o
main.o: main.c bellmanford.h
	$(CC) $(CFLAGS) -c main.c

# Rule for making bellmanford.o
bellmanford.o: bellmanford.c bellmanford.h
	$(CC) $(CFLAGS) -c bellmanford.c

# Clean the build directory
clean:
	rm -f $(TARGET) $(OBJS)
