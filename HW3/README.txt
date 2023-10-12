CSE6010 ASSIGNMENT 3 README
===========================

COMPILER AND OPERATING SYSTEM USED
---------------------------------
Compiler: gcc (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
Operating System: Ubuntu 20.04.4 LTS (WSL2 on Windows)
Kernel: Linux 5.10.102.1-microsoft-standard-WSL2 x86_64

DESCRIPTION
-----------
This assignment deals with identifying connected components in a graph, based on given edge labels. We implemented methods to read a graph from a file, perform depth-first search (DFS) on the graph, and count the number of connected components. The graph is represented using an adjacency list.

INSTRUCTIONS FOR COMPILING AND RUNNING
--------------------------------------

Using the Makefile (Recommended)
-------------------------------
1. Open a terminal and navigate to the folder containing the source code and Makefile (`~/desktop/6010/HW3`).
2. Run the command 'make' to compile the code.
3. Execute the compiled program by running './connections' followed by the filename and max label, e.g., `./connections data.txt 2`.

Manual Compilation for Main Program
-----------------------------------
1. Open a terminal and navigate to the folder containing the source code (`~/desktop/6010/HW3`).
2. Compile the `main.c` and `connections.c` files:
   ```bash
   gcc -Wall -Wextra -g -c main.c
   gcc -Wall -Wextra -g -c connections.c
