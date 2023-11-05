# The final executable name
particles: main.o particles.o
	gcc -o particles main.o particles.o -lm

# Compile main.c to main.o
# main.o depends on main.c and particles.h
main.o: main.c particles.h
	gcc -c -o main.o main.c

# Compile particles.c to particles.o
# particles.o depends on particles.c and particles.h
particles.o: particles.c particles.h
	gcc -c -o particles.o particles.c

# Clean target to remove object files and the executable
clean:
	rm -f *.o particles

