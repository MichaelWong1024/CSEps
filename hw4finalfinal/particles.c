#include "particles.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

Particle* particles;
Collision* collisions;
int numP;
double radius;
double lenX, lenY;
Collision lastCol;

// Macro definitions
#define minT(a, b) ((a) < (b) ? (a) : (b))
#define isEq(a, b, eps) (fabs((a) - (b)) < (eps))

/* Read the given file */
void readFromFile(char* fileName) {
    FILE* fp;
    fp = fopen(fileName, "r");
    
    fscanf(fp, "%d", &numP);
    fscanf(fp, "%lf", &radius);
    fscanf(fp, "%lf", &lenX);
    fscanf(fp, "%lf", &lenY);
    
    // allocate memory for the global variables
    particles = malloc(numP * sizeof(Particle)); // array of all particles
    collisions = malloc(numP * (numP + 1) / 2 * sizeof(Collision)); // array of all collisions
    
    // initialize the particle array
    for (int i = 0; i < numP; i++) {
        double x;
        double y;
        double vx;
        double vy;
        fscanf(fp, "%lf %lf %lf %lf", &x, &y, &vx, &vy);
        Particle p;
        p.id = i;
        p.x = x;
        p.y = y;
        p.vx = vx;
        p.vy = vy;
        p.tLast = 0;
        p.countP = 0;
        p.countW = 0;
        particles[i] = p;
    }
    
    fclose(fp);
}

/* Calculate the particle-wall collision based on the given particle */
Collision findColW(Particle p) {
    // A particle may collide on four sides of the box.
    // Here, we need to compare the t that the particle collides to each side,
    // and we choose the minimal t as the collision t.
    
    // initialize the t
    double tLeft = INFINITY;
    double tRight = INFINITY;
    double tTop = INFINITY;
    double tBottom = INFINITY;
    
    // calculate tLeft or tRight
    if (p.vx < 0) {
        // the particle moves to left
        tLeft = p.tLast + (radius - p.x) / p.vx;
    } else {
        // the particle moves to right
        tRight = p.tLast + ((lenX - radius) - p.x) / p.vx;
    }
    
    // calculate tTop or tBottom
    if (p.vy < 0) {
        // move to bottom
        tBottom = p.tLast + (radius - p.y) / p.vy;
    } else {
        // move to top
        tTop = p.tLast + ((lenY - radius) - p.y) / p.vy;
    }
    // choose the minimal t as the collision t
    double t = minT(tLeft, minT(tRight, minT(tTop, tBottom)));
    Collision collision;
    collision.p1 = p.id;
    collision.p2 = -1; // -1 represents the wall
    collision.t = t;
    return collision;
}

/* Calculate the particle-particle collision based on the two given particles */
Collision findColP(Particle p1, Particle p2, Collision lastCol) {
    // Since in each iteration, only part of the particles is updated,
    // the tLast of the two given particles might be different.
    // Here, we need to make the t of the two particles the same.
    if (!isEq(p1.tLast, p2.tLast, 0.0000001f)) {
        double tDiff = -1;
        if (p1.tLast > p2.tLast) {
            tDiff = p1.tLast - p2.tLast;
            p2.x += tDiff * p2.vx;
            p2.y += tDiff * p2.vy;
        } else {
            tDiff = p2.tLast - p1.tLast;
            p1.x += tDiff * p1.vx;
            p1.y += tDiff * p1.vy;
        }
    }
    
    // calculate the discriminant according to the formula in the guideline
    double t = INFINITY;
    double a = p1.x - p2.x;
    double b = p1.vx - p2.vx;
    double c = p1.y - p2.y;
    double d = p1.vy - p2.vy;
    double discriminant = pow(2 * a * b + 2 * c * d, 2) - 4 * (b * b + d * d) * (a * a + c * c - 4 * radius * radius);
    
    // based on the discriminent, calculate the collision t
    if (discriminant >= 0) {
        // calculate the smaller root as the collision t
        double res = (-2 * (a * b + c * d) - sqrt(discriminant)) / (2 * (b * b + d * d));
        // the root should be positive
        if (res >= 0) {
            t = res;
        }
    }
    Collision collision;
    collision.p1 = p1.id;
    collision.p2 = p2.id;
    collision.t = t + lastCol.t;
    
    return collision;
}

/* Initialize the collision array */
void initColArr(void) {
    // initialize the lastEarlist global variable
    Collision initialEarliest;
    initialEarliest.p1 = -1;
    initialEarliest.p2 = -1;
    initialEarliest.t = 0;
    lastCol = initialEarliest;
    
    // calculate particle-wall collisions
    // the first N items in the collision array are particle-wall collisions
    for (int i = 0; i < numP; i++) {
        Collision collision = findColW(particles[i]);
        collisions[i] = collision;
    }
    
    // calculate particle-particle collisions
    // the rest items in the collision array are particle-particle collisions
    int id = numP;
    for (int i = 0; i < numP; i++) {
        for (int j = i + 1; j < numP; j++) {
            collisions[id++] = findColP(particles[i], particles[j], lastCol);
        }
    }
}

/* sort the collision array using insertion sort */
void sortCollisions(void) {
    for (int i = 1; i < numP * (numP + 1) / 2; i++) {
        Collision temp = collisions[i];
        int j = i - 1;
        while (j >= 0 && collisions[j].t > temp.t) {
            collisions[j + 1] = collisions[j];
            j--;
        }
        collisions[j + 1] = temp;
    }
}

/* when particle-wall collision occurs, update the state of the collided particle */
void updateCollidedParticlePW(Particle particle, Collision earliest) {
    int id = particle.id;
    // calculate the location of the particle
    particles[id].x = particle.x + (earliest.t - particle.tLast) * particle.vx;
    particles[id].y = particle.y + (earliest.t - particle.tLast) * particle.vy;
    
    // if the particle collides on top or bottom wall
    if (isEq(particles[id].x, radius, 0.0000001f) || isEq(particles[id].x, lenX - radius, 0.0000001f)) {
        particles[id].vx = -particle.vx;
    }
    
    // if the particle collides on left or right wall
    if (isEq(particles[id].y, radius, 0.0000001f) || isEq(particles[id].y, lenY - radius, 0.0000001f)) {
        particles[id].vy = -particle.vy;
    }
    
    // update the last updated t
    particles[id].tLast = earliest.t;
    
    // update the particle-wall collision number of the particle
    particles[id].countW++;
}

/* when particle-particle collision occurs, update the state of the collided particles */
void updateCollidedParticlePP(Particle p1, Particle p2, Collision earliest) {
    int id1 = p1.id;
    int id2 = p2.id;
    
    // update the location of the p1
    particles[id1].x = p1.x + (earliest.t - p1.tLast) * p1.vx;
    particles[id1].y = p1.y + (earliest.t - p1.tLast) * p1.vy;
    
    // update the location of the p2
    particles[id2].x = p2.x + (earliest.t - p2.tLast) * p2.vx;
    particles[id2].y = p2.y + (earliest.t - p2.tLast) * p2.vy;
    
    // swap the x-axis velocity
    double temp = particles[id1].vx;
    particles[id1].vx = particles[id2].vx;
    particles[id2].vx = temp;
    
    // swap the y-axis velocity
    temp = particles[id1].vy;
    particles[id1].vy = particles[id2].vy;
    particles[id2].vy = temp;
    
    // update the last updated t
    particles[id1].tLast = earliest.t;
    particles[id2].tLast = earliest.t;
    
    // update the particle-particle collision number of the particles
    particles[id1].countP++;
    particles[id2].countP++;
}

/* Process the earliest collision in the collision array */
void processEarliest(void) {
    Collision earliest = collisions[0];
    // record the earliest collision before the collision array is updated
    lastCol = earliest;
    if (earliest.p2 == -1) {
        // particle-wall collision
        Particle particle = particles[earliest.p1];
        updateCollidedParticlePW(particle, earliest);
    } else {
        // particle-particle collision
        Particle p1 = particles[earliest.p1];
        Particle p2 = particles[earliest.p2];
        updateCollidedParticlePP(p1, p2, earliest);
    }
}

/* After processing the earliest collision, find all the collision involving the particles in the earliest collision */
void updateCollisions(void) {
    Collision earliest = lastCol;
    if (earliest.p2 == -1) {
        // the earliest collision is particle-wall collision
        Particle particle = particles[earliest.p1];
        for (int i = 0; i < numP * (numP + 1) / 2; i++) {
            if (collisions[i].p1 == particle.id && collisions[i].p2 == -1) {
                // the scanned collision is a particle-wall collision
                Collision newCollision = findColW(particle);
                collisions[i].t = newCollision.t;
            } else if (collisions[i].p1 == particle.id && collisions[i].p2 != -1) {
                // the scanned collision is a particle-particle collision
                Collision newCollision = findColP(particle, particles[collisions[i].p2], lastCol);
                collisions[i].t = newCollision.t;
            } else if (collisions[i].p1 != particle.id && collisions[i].p2 == particle.id) {
                // the scanned collision is a particle-wall collision
                Collision newCollision = findColP(particles[collisions[i].p1], particle, lastCol);
                collisions[i].t = newCollision.t;
            }
        }
    } else {
        // the earliest collision is particle-particle collision
        Particle p1 = particles[earliest.p1];
        Particle p2 = particles[earliest.p2];
        for (int i = 0; i < numP * (numP + 1) / 2; i++) {
            if (collisions[i].p2 != -1) {
                // the scanned collision is a particle-particle collision
                if (collisions[i].p1 == p1.id || collisions[i].p1 == p2.id || collisions[i].p2 == p1.id || collisions[i].p2 == p2.id) {
                    Collision newCollision = findColP(particles[collisions[i].p1], particles[collisions[i].p2], lastCol);
                    collisions[i].t = newCollision.t;
                }
            } else {
                // the scanned collision is a particle-wall collision
                if (collisions[i].p1 == p1.id || collisions[i].p1 == p2.id) {
                    Collision newCollision = findColW(particles[collisions[i].p1]);
                    collisions[i].t = newCollision.t;
                }
            }
        }
    }
}

/* After the last collision is processed, update the positions of all the particles to the last collision t */
void updateAllParticles(void) {
    int p1 = lastCol.p1;
    int p2 = lastCol.p2;
    double t = lastCol.t;
    for (int id = 0; id < numP; id++) {
        if (id != p1 && id != p2) {
            particles[id].x = particles[id].x + (t - particles[id].tLast) * particles[id].vx;
            particles[id].y = particles[id].y + (t - particles[id].tLast) * particles[id].vy;
            particles[id].tLast = t;
        }
    }
}

/* free the allocated memory */
void finalize(void) {
    free(particles);
    free(collisions);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("Wrong number of arguments!\n");
        return 0;
    }
    char* filename = argv[1];
    double endTime = atof(argv[2]);
    // Initialize the particle array and the collision array
    readFromFile(filename);
    initColArr();
    sortCollisions();
    // Iterate until the earliest collision t larger than the t given by the user
    while (collisions[0].t <= endTime) {
        processEarliest();
        updateCollisions();
        sortCollisions();
    }
    // After the last collision is processed, update the positions of all the particles to the last collision t
    updateAllParticles();
    // Update the positions of all the particles to the end t given by the user
    for (int id = 0; id < numP; id++) {
        particles[id].x = particles[id].x + (endTime - particles[id].tLast) * particles[id].vx;
        particles[id].y = particles[id].y + (endTime - particles[id].tLast) * particles[id].vy;
        printf("%.6f, %.6f, %d, %d\n", particles[id].x, particles[id].y, particles[id].countW, particles[id].countP);
    }
    finalize();
    return 1;
}
