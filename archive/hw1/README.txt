CSE6010 HOMEWORK 1 README
=========================

COMPILER AND OPERATING SYSTEM USED
---------------------------------
Compiler: gcc (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
Operating System: Ubuntu 20.04.4 LTS (WSL2 on Windows)
Kernel: Linux 5.10.102.1-microsoft-standard-WSL2 x86_64

INSTRUCTIONS FOR COMPILING AND RUNNING
--------------------------------------

Using the Makefile (Recommended)
--------------------------------
1. Open a terminal and navigate to the folder containing the source code and Makefile (`~/desktop/CSE6010/HW1`).
2. Run the command 'make' to compile the code.
3. Execute the output file by running './solution'.

Manual Compilation for Main Program
-----------------------------------
1. Open a terminal and navigate to the folder containing the source code.
2. Compile each .c file into an object file using the following commands:
    gcc -c main.c -o main.o -Wall -Wsign-compare -Wwrite-strings -Wtype-limits -Werror -fsanitize=address
    gcc -c solution.c -o solution.o -Wall -Wsign-compare -Wwrite-strings -Wtype-limits -Werror -fsanitize=address
3. Link the object files to create the executable:
    gcc main.o solution.o -o solution -lm -Wall -Wsign-compare -Wwrite-strings -Wtype-limits -Werror -fsanitize=address
4. Execute the output file by running './solution'.

Manual Compilation for Stress Test
----------------------------------
1. Open a terminal and navigate to the folder containing the stress test source code (`stress_test.c`) and the main solution file (`solution.c`).
2. Compile the stress test and main solution into an executable using the following command:
    gcc -o stress_test solution.c stress_test.c -lm -Wall -Wsign-compare -Wwrite-strings -Wtype-limits -Werror -fsanitize=address
3. Execute the stress test by running './stress_test'.
