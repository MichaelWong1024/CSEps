CSE6010 ASSIGNMENT 2 README
===========================

COMPILER AND OPERATING SYSTEM USED
---------------------------------
Compiler: gcc (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
Operating System: Ubuntu 20.04.4 LTS (WSL2 on Windows)
Kernel: Linux 5.10.102.1-microsoft-standard-WSL2 x86_64

INSTRUCTIONS FOR COMPILING AND RUNNING
--------------------------------------

Using the Makefile (Recommended)
-------------------------------
1. Open a terminal and navigate to the folder containing the source code and Makefile (`~/desktop/6010/HW2`).
2. Run the command 'make' to compile the code.
3. Execute the compiled program by running './pq'.

Manual Compilation for Main Program
-----------------------------------
1. Open a terminal and navigate to the folder containing the source code (`~/desktop/6010/HW2`).
2. Compile the `main.c` and `pq.c` files:
   gcc -ggdb3 -c main.c
   gcc -ggdb3 -c pq.c
3. Link the object files to create the executable:
   gcc -ggdb3 -o pq main.o pq.o -lm
4. Execute the output file by running './pq'.
