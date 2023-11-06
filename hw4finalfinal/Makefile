# The final executable will be named 'particle'
particle: main.o particles.o
	gcc -o particles main.o particles.o -lm

# This will compile any .c to .o
%.o: %.c particles.h
	gcc -c -o $@ $<

# This is the clean command to remove all object files and the executable
clean:
	rm -f *.o particles

