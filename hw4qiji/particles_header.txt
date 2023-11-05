#ifndef particles_h
#define particles_h

#include <stdio.h>
typedef struct particle {
    int index;
    double currX;
    double currY;
    double currVx;
    double currVy;
    double lastTime;
    int numCollisionPW;
    int numCollisionPP;
} Particle;

typedef struct collision {
    int particle1;
    int particle2;
    double time;
} Collision;

void readFromFile(char* fileName);
double min(double x, double y);
int areDoublesEqual(double x, double y, double epsilon);
Collision calCollisionPW(Particle p);
Collision calCollisionPP(Particle p1, Particle p2, Collision lastEarliest);
void calCollisions(void);
void sortCollisions(void);
void calInitialCollisions(void);
void sortCollisions(void);
void updateCollidedParticlePW(Particle particle, Collision earliest);
void updateCollidedParticlePP(Particle particle1, Particle particle2, Collision earliest);
void processEarliest(void);
void updateCollisions(void);
void updateAllParticles(void);
void finalize(void);
#endif /* particles_h */
