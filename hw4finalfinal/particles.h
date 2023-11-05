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

typedef struct collision {
    int particle1;
    int particle2;
    double time;
} Collision;

// typedef struct {
//     int p1, p2;       // Colliding particles' IDs
//     double t;         // Collision time
// } Collision;

void readFromFile(char* fileName);

double min(double x, double y);
int areDoublesEqual(double x, double y, double epsilon);

Collision calCollisionPW(Particle p);
Collision calCollisionPP(Particle p1, Particle p2, Collision lastEarliest);

void calCollisions(void);
void sortCollisions(void);
void calInitialCollisions(void);
void updateCollidedParticlePW(Particle particle, Collision earliest);
void updateCollidedParticlePP(Particle particle1, Particle particle2, Collision earliest);
void processEarliest(void);
void updateCollisions(void);
void updateAllParticles(void);
void finalize(void);

#endif /* PARTICLES_H */
