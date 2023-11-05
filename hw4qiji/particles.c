#include "particles.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

Particle* particles;
Collision* collisions;
int N;
double radius;
double lenX;
double lenY;
Collision lastEarliest;

/* Read the given file */
void readFromFile(char* fileName) {
    FILE* fp;
    fp = fopen(fileName, "r");
    
    fscanf(fp, "%d", &N);
    fscanf(fp, "%lf", &radius);
    fscanf(fp, "%lf", &lenX);
    fscanf(fp, "%lf", &lenY);
    
    // allocate memory for the global variables
    particles = malloc(N * sizeof(Particle)); // array of all particles
    collisions = malloc(N * (N + 1) / 2 * sizeof(Collision)); // array of all collisions
    
    // initialize the particle array
    for (int i = 0; i < N; i++) {
        double currX;
        double currY;
        double currVx;
        double currVy;
        fscanf(fp, "%lf %lf %lf %lf", &currX, &currY, &currVx, &currVy);
        Particle p;
        p.index = i;
        p.currX = currX;
        p.currY = currY;
        p.currVx = currVx;
        p.currVy = currVy;
        p.lastTime = 0;
        p.numCollisionPP = 0;
        p.numCollisionPW = 0;
        particles[i] = p;
    }
    
    fclose(fp);
}

/* Helper function to compare two doubles */
double min(double x, double y) {
    if (x < y) {
        return x;
    } else {
        return y;
    }
}

/* Helper function to check whether two doubles equals within certain precision epsilon */
int areDoublesEqual(double x, double y, double epsilon) {
    return fabs(x - y) < epsilon;
}

/* Calculate the particle-wall collision based on the given particle */
Collision calCollisionPW(Particle p) {
    // A particle may collide on four sides of the box.
    // Here, we need to compare the time that the particle collides to each side,
    // and we choose the minimal time as the collision time.
    
    // initialize the time
    double tLeft = INFINITY;
    double tRight = INFINITY;
    double tTop = INFINITY;
    double tBottom = INFINITY;
    
    // calculate tLeft or tRight
    if (p.currVx < 0) {
        // the particle moves to left
        tLeft = p.lastTime + (radius - p.currX) / p.currVx;
    } else {
        // the particle moves to right
        tRight = p.lastTime + ((lenX - radius) - p.currX) / p.currVx;
    }
    
    // calculate tTop or tBottom
    if (p.currVy < 0) {
        // move to bottom
        tBottom = p.lastTime + (radius - p.currY) / p.currVy;
    } else {
        // move to top
        tTop = p.lastTime + ((lenY - radius) - p.currY) / p.currVy;
    }
    // choose the minimal time as the collision time
    double time = min(tLeft, min(tRight, min(tTop, tBottom)));
    Collision collision;
    collision.particle1 = p.index;
    collision.particle2 = -1; // -1 represents the wall
    collision.time = time;
    return collision;
}

/* Calculate the particle-particle collision based on the two given particles */
Collision calCollisionPP(Particle p1, Particle p2, Collision lastEarliest) {
    // Since in each iteration, only part of the particles is updated,
    // the lastTime of the two given particles might be different.
    // Here, we need to make the time of the two particles the same.
    if (!areDoublesEqual(p1.lastTime, p2.lastTime, 0.0000001f)) {
        double timeDiff = -1;
        if (p1.lastTime > p2.lastTime) {
            timeDiff = p1.lastTime - p2.lastTime;
            p2.currX += timeDiff * p2.currVx;
            p2.currY += timeDiff * p2.currVy;
        } else {
            timeDiff = p2.lastTime - p1.lastTime;
            p1.currX += timeDiff * p1.currVx;
            p1.currY += timeDiff * p1.currVy;
        }
    }
    
    // calculate the discriminant according to the formula in the guideline
    double time = INFINITY;
    double a = p1.currX - p2.currX;
    double b = p1.currVx - p2.currVx;
    double c = p1.currY - p2.currY;
    double d = p1.currVy - p2.currVy;
    double discriminant = pow(2 * a * b + 2 * c * d, 2) - 4 * (b * b + d * d) * (a * a + c * c - 4 * radius * radius);
    
    // based on the discriminent, calculate the collision time
    if (discriminant >= 0) {
        // calculate the smaller root as the collision time
        double res = (-2 * (a * b + c * d) - sqrt(discriminant)) / (2 * (b * b + d * d));
        // the root should be positive
        if (res >= 0) {
            time = res;
        }
    }
    Collision collision;
    collision.particle1 = p1.index;
    collision.particle2 = p2.index;
    collision.time = time + lastEarliest.time;
    
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
