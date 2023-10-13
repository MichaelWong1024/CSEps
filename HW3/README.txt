CSE6010 ASSIGNMENT 3 README
===========================

COMPILER AND OPERATING SYSTEM USED
---------------------------------
Compiler: gcc (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
Operating System: Ubuntu 20.04.4 LTS (WSL2 on Windows)
Kernel: Linux 5.10.102.1-microsoft-standard-WSL2 x86_64

DESCRIPTION
-----------
This program focuses on the identification of connected components in a graph, with the additional complexity of varying edge weights.
The program reads a graph from a file, performs depth-first search (DFS) to determine the connected components, and presents the count of these components.
The graph is represented using an adjacency list.

TESTING
-------
Stress Testing Approach:
- Comprehensive testing of complex graph structures with varying connected components.
- Ensured detection of distinct connected components with alternating edge weights.
- Validation of correct weight interpretation.

Validation Criteria:
- Depth-First Search must correctly identify separate connected components in the graph.
- Accurate reading and processing of edge weights.
- Handling of repeated edges and unique vertex identifiers.

Results:
All tests, including the advanced graph structure test (`connections_advanced_test.txt`), were successfully executed. The program efficiently determines connected components and respects the given edge weights in diverse graph structures.

INSTRUCTIONS FOR COMPILING AND RUNNING
--------------------------------------

Using the Makefile (Recommended)
-------------------------------
1. Open a terminal and navigate to the folder containing the source code and Makefile (`~/desktop/6010/HW3`).
2. Run the command 'make' to compile the code.
3. Execute the compiled program by running './connections' followed by the filename and max label, e.g., `./connections connections_advanced_test.txt 2`.

Manual Compilation for Main Program
-----------------------------------
1. Open a terminal and navigate to the folder containing the source code (`~/desktop/6010/HW3`).
2. Compile the `main.c` and `connections.c` files:
   ```bash
   gcc -Wall -Wextra -g -c main.c
   gcc -Wall -Wextra -g -c connections.c
