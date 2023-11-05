#include "particles.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Global Variables
Particle *particles;
Collision *collisions;
int nParticles;     // Number of particles
double radius;      // Radius of particles
double lenX, lenY;  // Dimensions of the container
Collision lastEarliest;

// Macro definitions
#define minT(a, b) ((a) < (b) ? (a) : (b))
#define isEq(a, b, eps) (fabs((a) - (b)) < (eps))

// Reads initial configuration from file
void read(const char* file) {
    FILE* fp = fopen(file, "r");
    if (fp == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }
    
    // Read parameters and allocate memory
    fscanf(fp, "%d %lf %lf %lf", &nParticles, &radius, &lenX, &lenY);
    particles = malloc(nParticles * sizeof(Particle));
    collisions = malloc(nParticles * (nParticles + 1) / 2 * sizeof(Collision));

    // Read particle initial conditions
    for (int i = 0; i < nParticles; ++i) {
        Particle p = {i, 0, 0, 0, 0, 0, 0, 0};
        fscanf(fp, "%lf %lf %lf %lf", &p.x, &p.y, &p.vx, &p.vy);
        particles[i] = p;
    }
    
    fclose(fp);
}

/* Calculate the particle-wall collision based on the given particle */
Collision findWC(Particle p) {
    // A particle may collide on four sides of the box.
    // Here, we need to compare the time at which the particle collides with each side,
    // and we choose the minimal time as the collision time.
    
    // initialize the time
    double tLeft = INFINITY;
    double tRight = INFINITY;
    double tTop = INFINITY;
    double tBottom = INFINITY;
    
    // calculate tLeft or tRight
    if (p.vx < 0) {
        // the particle moves to the left
        tLeft = (radius - p.x) / p.vx;
    } else if (p.vx > 0) {
        // the particle moves to the right
        tRight = ((lenX - radius) - p.x) / p.vx;
    }
    
    // calculate tTop or tBottom
    if (p.vy < 0) {
        // the particle moves to the bottom
        tBottom = (radius - p.y) / p.vy;
    } else if (p.vy > 0) {
        // the particle moves to the top
        tTop = ((lenY - radius) - p.y) / p.vy;
    }

    // Choose the minimal time as the collision time, relative to the last update time
    double time = p.tLast + minT(minT(tLeft, tRight), minT(tTop, tBottom));

    Collision collision;
    collision.p1 = p.id;   // Updated from 'particle1' to 'p1' to match the Collision structure
    collision.p2 = -1;     // -1 represents the wall
    collision.t = time;    // Updated from 'time' to 't' to match the Collision structure
    return collision;
}

/* Calculate the particle-particle collision based on the two given particles */
Collision calCollisionPP(Particle p1, Particle p2, Collision lastEarliest) {
    // Since in each iteration, only part of the particles is updated,
    // the tLast of the two given particles might be different.
    // Here, we need to make the time of the two particles the same.
    if (!areDoublesEqual(p1.tLast, p2.tLast, 0.0000001f)) {
        double timeDiff = -1;
        if (p1.tLast > p2.tLast) {
            timeDiff = p1.tLast - p2.tLast;
            p2.x += timeDiff * p2.vx;
            p2.y += timeDiff * p2.vy;
        } else {
            timeDiff = p2.tLast - p1.tLast;
            p1.x += timeDiff * p1.vx;
            p1.y += timeDiff * p1.vy;
        }
    }
    
    // calculate the discriminant according to the formula in the guideline
    double time = INFINITY;
    double dx = p1.x - p2.x;
    double dvx = p1.vx - p2.vx;
    double dy = p1.y - p2.y;
    double dvy = p1.vy - p2.vy;
    double discriminant = pow(2 * dx * dvx + 2 * dy * dvy, 2) - 4 * (dvx * dvx + dvy * dvy) * (dx * dx + dy * dy - 4 * radius * radius);
    
    // based on the discriminant, calculate the collision time
    if (discriminant >= 0) {
        // calculate the smaller root as the collision time
        double root = (-2 * (dx * dvx + dy * dvy) - sqrt(discriminant)) / (2 * (dvx * dvx + dvy * dvy));
        // the root should be positive
        if (root >= 0) {
            time = root;
        }
    }
    Collision collision;
    collision.p1 = p1.id;  // corrected to id
    collision.p2 = p2.id;  // corrected to id
    collision.t = time + lastEarliest.t;  // corrected to t
    
    return collision;
}

/* Initialize the collision array */
void calInitialCollisions(void) {
    // initialize the lastEarlist global variable
    Collision initialEarliest;
    initialEarliest.particle1 = -1;
    initialEarliest.particle2 = -1;
    initialEarliest.time = 0;
    lastEarliest = initialEarliest;
    
    // calculate particle-wall collisions
    // the first N items in the collision array are particle-wall collisions
    for (int i = 0; i < N; i++) {
        Collision collision = calCollisionPW(particles[i]);
        collisions[i] = collision;
    }
    
    // calculate particle-particle collisions
    // the rest items in the collision array are particle-particle collisions
    int index = N;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            collisions[index++] = calCollisionPP(particles[i], particles[j], lastEarliest);
        }
    }
}

/* sort the collision array using insertion sort */
void sortCollisions(void) {
    for (int i = 1; i < N * (N + 1) / 2; i++) {
        Collision temp = collisions[i];
        int j = i - 1;
        while (j >= 0 && collisions[j].time > temp.time) {
            collisions[j + 1] = collisions[j];
            j--;
        }
        collisions[j + 1] = temp;
    }
}

/* when particle-wall collision occurs, update the state of the collided particle */
void updateCollidedParticlePW(Particle particle, Collision earliest) {
    int index = particle.index;
    // calculate the location of the particle
    particles[index].currX = particle.currX + (earliest.time - particle.lastTime) * particle.currVx;
    particles[index].currY = particle.currY + (earliest.time - particle.lastTime) * particle.currVy;
    
    // if the particle collides on top or bottom wall
    if (areDoublesEqual(particles[index].currX, radius, 0.0000001f) || areDoublesEqual(particles[index].currX, lenX - radius, 0.0000001f)) {
        particles[index].currVx = -particle.currVx;
    }
    
    // if the particle collides on left or right wall
    if (areDoublesEqual(particles[index].currY, radius, 0.0000001f) || areDoublesEqual(particles[index].currY, lenY - radius, 0.0000001f)) {
        particles[index].currVy = -particle.currVy;
    }
    
    // update the last updated time
    particles[index].lastTime = earliest.time;
    
    // update the particle-wall collision number of the particle
    particles[index].numCollisionPW++;
}

/* when particle-particle collision occurs, update the state of the collided particles */
void updateCollidedParticlePP(Particle particle1, Particle particle2, Collision earliest) {
    int index1 = particle1.index;
    int index2 = particle2.index;
    
    // update the location of the particle1
    particles[index1].currX = particle1.currX + (earliest.time - particle1.lastTime) * particle1.currVx;
    particles[index1].currY = particle1.currY + (earliest.time - particle1.lastTime) * particle1.currVy;
    
    // update the location of the particle2
    particles[index2].currX = particle2.currX + (earliest.time - particle2.lastTime) * particle2.currVx;
    particles[index2].currY = particle2.currY + (earliest.time - particle2.lastTime) * particle2.currVy;
    
    // swap the x-axis velocity
    double temp = particles[index1].currVx;
    particles[index1].currVx = particles[index2].currVx;
    particles[index2].currVx = temp;
    
    // swap the y-axis velocity
    temp = particles[index1].currVy;
    particles[index1].currVy = particles[index2].currVy;
    particles[index2].currVy = temp;
    
    // update the last updated time
    particles[index1].lastTime = earliest.time;
    particles[index2].lastTime = earliest.time;
    
    // update the particle-particle collision number of the particles
    particles[index1].numCollisionPP++;
    particles[index2].numCollisionPP++;
}

/* Process the earliest collision in the collision array */
void processEarliest(void) {
    Collision earliest = collisions[0];
    // record the earliest collision before the collision array is updated
    lastEarliest = earliest;
    if (earliest.particle2 == -1) {
        // particle-wall collision
        Particle particle = particles[earliest.particle1];
        updateCollidedParticlePW(particle, earliest);
    } else {
        // particle-particle collision
        Particle particle1 = particles[earliest.particle1];
        Particle particle2 = particles[earliest.particle2];
        updateCollidedParticlePP(particle1, particle2, earliest);
    }
}

/* After processing the earliest collision, find all the collision involving the particles in the earliest collision */
void updateCollisions(void) {
    Collision earliest = lastEarliest;
    if (earliest.particle2 == -1) {
        // the earliest collision is particle-wall collision
        Particle particle = particles[earliest.particle1];
        for (int i = 0; i < N * (N + 1) / 2; i++) {
            if (collisions[i].particle1 == particle.index && collisions[i].particle2 == -1) {
                // the scanned collision is a particle-wall collision
                Collision newCollision = calCollisionPW(particle);
                collisions[i].time = newCollision.time;
            } else if (collisions[i].particle1 == particle.index && collisions[i].particle2 != -1) {
                // the scanned collision is a particle-particle collision
                Collision newCollision = calCollisionPP(particle, particles[collisions[i].particle2], lastEarliest);
                collisions[i].time = newCollision.time;
            } else if (collisions[i].particle1 != particle.index && collisions[i].particle2 == particle.index) {
                // the scanned collision is a particle-wall collision
                Collision newCollision = calCollisionPP(particles[collisions[i].particle1], particle, lastEarliest);
                collisions[i].time = newCollision.time;
            }
        }
    } else {
        // the earliest collision is particle-particle collision
        Particle particle1 = particles[earliest.particle1];
        Particle particle2 = particles[earliest.particle2];
        for (int i = 0; i < N * (N + 1) / 2; i++) {
            if (collisions[i].particle2 != -1) {
                // the scanned collision is a particle-particle collision
                if (collisions[i].particle1 == particle1.index || collisions[i].particle1 == particle2.index || collisions[i].particle2 == particle1.index || collisions[i].particle2 == particle2.index) {
                    Collision newCollision = calCollisionPP(particles[collisions[i].particle1], particles[collisions[i].particle2], lastEarliest);
                    collisions[i].time = newCollision.time;
                }
            } else {
                // the scanned collision is a particle-wall collision
                if (collisions[i].particle1 == particle1.index || collisions[i].particle1 == particle2.index) {
                    Collision newCollision = calCollisionPW(particles[collisions[i].particle1]);
                    collisions[i].time = newCollision.time;
                }
            }
        }
    }
}

/* After the last collision is processed, update the positions of all the particles to the last collision time */
void updateAllParticles(void) {
    int particle1 = lastEarliest.particle1;
    int particle2 = lastEarliest.particle2;
    double time = lastEarliest.time;
    for (int index = 0; index < N; index++) {
        if (index != particle1 && index != particle2) {
            particles[index].currX = particles[index].currX + (time - particles[index].lastTime) * particles[index].currVx;
            particles[index].currY = particles[index].currY + (time - particles[index].lastTime) * particles[index].currVy;
            particles[index].lastTime = time;
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
    calInitialCollisions();
    sortCollisions();
    // Iterate until the earliest collision time larger than the time given by the user
    while (collisions[0].time <= endTime) {
        processEarliest();
        updateCollisions();
        sortCollisions();
    }
    // After the last collision is processed, update the positions of all the particles to the last collision time
    updateAllParticles();
    // Update the positions of all the particles to the end time given by the user
    for (int index = 0; index < N; index++) {
        particles[index].currX = particles[index].currX + (endTime - particles[index].lastTime) * particles[index].currVx;
        particles[index].currY = particles[index].currY + (endTime - particles[index].lastTime) * particles[index].currVy;
        printf("%.6f, %.6f, %d, %d\n", particles[index].currX, particles[index].currY, particles[index].numCollisionPW, particles[index].numCollisionPP);
    }
    finalize();
    return 1;
}
