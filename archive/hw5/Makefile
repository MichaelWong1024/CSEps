CC = gcc
CFLAGS = -Wall -g
LDFLAGS = -fopenmp

TARGET = bellmanford

OBJS = main.o bellmanford.o

all: $(TARGET)

$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

main.o: main.c bellmanford.h
	$(CC) $(CFLAGS) -c main.c

bellmanford.o: bellmanford.c bellmanford.h
	$(CC) $(CFLAGS) -c bellmanford.c

clean:
	rm -f $(TARGET) $(OBJS)

.PHONY: all clean

