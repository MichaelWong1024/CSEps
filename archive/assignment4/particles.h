#ifndef PARTICLES_H
#define PARTICLES_H

#include <stdio.h>

// Particle and Collision Structures
typedef struct {
    int id;
    double x, y;      // Position
    double vx, vy;    // Velocity
    double tLast;     // Last update time
    int countW, countP;      // Wall, Particle collisions count
} Particle;

typedef struct {
    int p1, p2;       // Colliding particles' IDs
    double t;         // Collision time
} Collision;

extern Particle* particles;
extern Collision* collisions;
extern int numP;
extern double radius;
extern double lenX;
extern double lenY;
extern Collision lastCol;

void read(const char* fileName);

// Collision calculation functions
Collision findColW(Particle p); // Find wall collision
Collision findColP(Particle p1, Particle p2, Collision lastEarliest); // Find particle collision

// Collision handling functions
void initColArr(); // Compute all collisions into an array
void sortColArr(); // Sort collisions array
void updatePW(Particle particle, Collision earliest); // Update wall collision
void updatePP(Particle particle1, Particle particle2, Collision earliest); // Update particle collision
void jumpCol(); // Process the next collision
void updateColArr(); // Update collision list
void updateParArr(); // Update all particles state
void finalize(); // Finalize the simulation

#endif /* PARTICLES_H */
