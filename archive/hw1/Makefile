solution: main.o solution.o
	gcc -o solution main.o solution.o -lm

%.o: %.c
	gcc -c -o $@ $<

clean:
	rm -f *.o solution
