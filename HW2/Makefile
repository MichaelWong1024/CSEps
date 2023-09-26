FLAGS = -ggdb3

solution: main.o pq.o
	        gcc $(FLAGS) -o pq main.o pq.o -lm

%.o: %.c
	        gcc $(FLAGS) -c -o $@ $<

clean:
	        rm -f *.o pq
