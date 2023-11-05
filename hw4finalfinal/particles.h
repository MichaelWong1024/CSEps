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

void readFromFile(char* fileName);

// Collision calculation functions
Collision findColW(Particle p);
Collision findColP(Particle p1, Particle p2, Collision lastEarliest);

void initColArr(void);

void calCollisions(void);
void sortCollisions(void);
void updateCollidedParticlePW(Particle particle, Collision earliest);
void updateCollidedParticlePP(Particle particle1, Particle particle2, Collision earliest);
void processEarliest(void);
void updateCollisions(void);
void updateAllParticles(void);
void finalize(void);

#endif /* PARTICLES_H */
