#ifndef PARTICLES_H
#define PARTICLES_H

#include <stdio.h>

// Particle and Collision Structures
typedef struct {
    int id;
    double x, y;      // Position
    double vx, vy;    // Velocity
    double tLast;     // Last update time
    int colWall;      // Wall collisions count
    int colPart;      // Particle collisions count
} Particle;

typedef struct {
    int p1, p2;       // Colliding particles' IDs
    double t;         // Collision time
} Collision;

void read(char* file);               // Read data from file

// Collision calculation functions
Collision findWC(Particle p);      // Find wall collision
Collision findPC(Particle p1, Particle p2, Collision last); // Find particle collision

// Collision handling functions
void computeAll();                        // Compute all collisions
void sortCols();                          // Sort collisions
void initializeC();                          // Initialize the collision array
void updateWC(Particle* p, Collision* col); // Update wall collision
void updatePC(Particle* p1, Particle* p2, const Collision* col); // Update particle collision
void processNextC();                       // Process the next collision
void updateColsL();                       // Refresh collision list
void updateAllP();                         // Update all particles
void finalize();                            // Finalize the simulation

#endif /* PARTICLES_H */
