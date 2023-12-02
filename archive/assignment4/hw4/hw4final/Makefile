solution: particles.o
	gcc -o particles particles.o -lm

%.o: %.c
	gcc -c -o $@ $<

clean:
	rm -f *.o particles
