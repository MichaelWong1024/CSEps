#include <stdio.h>
#include <stdlib.h>
#include "particles.h"

/* Main simulation control flow */
int main(int argc, char* argv[]) {
    // Check for correct command-line argument count
    if (argc != 3) {
        printf("Usage: %s <filename> <endTime>\n", argv[0]);
        return 0; // Return 0 for unsuccessful execution
    }

    // Initialization from input
    double endTime = atof(argv[2]);
    read(argv[1]); // Read input file
    initColArr();  // Initialize collision array
    sortColArr();  // Sort collision array

    // Main simulation loop
    while (collisions[0].t <= endTime) {
        jumpCol();      // Process earliest collision
        updateColArr(); // Update collision array
        sortColArr();   // Sort collision array after update
    }

    // Update particles' positions to the last collision time
    updateParArr();

    // Update particles' positions to end time and print results
    for (int id = 0; id < numP; id++) {
        double deltaX = (endTime - particles[id].tLast) * particles[id].vx;
        double deltaY = (endTime - particles[id].tLast) * particles[id].vy;
        particles[id].x += deltaX;
        particles[id].y += deltaY;
        printf("%.6f, %.6f, %d, %d\n", 
            particles[id].x, particles[id].y, 
            particles[id].countW, particles[id].countP);
    }

    finalize(); // Clean-up resources
    return 1; // Return 1 for successful execution
}

