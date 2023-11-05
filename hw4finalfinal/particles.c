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
    // choose the minimal time as the collision time
    double time = min(tLeft, min(tRight, min(tTop, tBottom)));
    Collision collision;
    collision.particle1 = p.id;
    collision.particle2 = -1; // -1 represents the wall
    collision.time = time;
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
    double a = p1.x - p2.x;
    double b = p1.vx - p2.vx;
    double c = p1.y - p2.y;
    double d = p1.vy - p2.vy;
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
    collision.particle1 = p1.id;
    collision.particle2 = p2.id;
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
    int id = N;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            collisions[id++] = calCollisionPP(particles[i], particles[j], lastEarliest);
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
    int id = particle.id;
    // calculate the location of the particle
    particles[id].x = particle.x + (earliest.time - particle.tLast) * particle.vx;
    particles[id].y = particle.y + (earliest.time - particle.tLast) * particle.vy;
    
    // if the particle collides on top or bottom wall
    if (areDoublesEqual(particles[id].x, radius, 0.0000001f) || areDoublesEqual(particles[id].x, lenX - radius, 0.0000001f)) {
        particles[id].vx = -particle.vx;
    }
    
    // if the particle collides on left or right wall
    if (areDoublesEqual(particles[id].y, radius, 0.0000001f) || areDoublesEqual(particles[id].y, lenY - radius, 0.0000001f)) {
        particles[id].vy = -particle.vy;
    }
    
    // update the last updated time
    particles[id].tLast = earliest.time;
    
    // update the particle-wall collision number of the particle
    particles[id].countW++;
}

/* when particle-particle collision occurs, update the state of the collided particles */
void updateCollidedParticlePP(Particle particle1, Particle particle2, Collision earliest) {
    int id1 = particle1.id;
    int id2 = particle2.id;
    
    // update the location of the particle1
    particles[id1].x = particle1.x + (earliest.time - particle1.tLast) * particle1.vx;
    particles[id1].y = particle1.y + (earliest.time - particle1.tLast) * particle1.vy;
    
    // update the location of the particle2
    particles[id2].x = particle2.x + (earliest.time - particle2.tLast) * particle2.vx;
    particles[id2].y = particle2.y + (earliest.time - particle2.tLast) * particle2.vy;
    
    // swap the x-axis velocity
    double temp = particles[id1].vx;
    particles[id1].vx = particles[id2].vx;
    particles[id2].vx = temp;
    
    // swap the y-axis velocity
    temp = particles[id1].vy;
    particles[id1].vy = particles[id2].vy;
    particles[id2].vy = temp;
    
    // update the last updated time
    particles[id1].tLast = earliest.time;
    particles[id2].tLast = earliest.time;
    
    // update the particle-particle collision number of the particles
    particles[id1].countP++;
    particles[id2].countP++;
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
            if (collisions[i].particle1 == particle.id && collisions[i].particle2 == -1) {
                // the scanned collision is a particle-wall collision
                Collision newCollision = calCollisionPW(particle);
                collisions[i].time = newCollision.time;
            } else if (collisions[i].particle1 == particle.id && collisions[i].particle2 != -1) {
                // the scanned collision is a particle-particle collision
                Collision newCollision = calCollisionPP(particle, particles[collisions[i].particle2], lastEarliest);
                collisions[i].time = newCollision.time;
            } else if (collisions[i].particle1 != particle.id && collisions[i].particle2 == particle.id) {
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
                if (collisions[i].particle1 == particle1.id || collisions[i].particle1 == particle2.id || collisions[i].particle2 == particle1.id || collisions[i].particle2 == particle2.id) {
                    Collision newCollision = calCollisionPP(particles[collisions[i].particle1], particles[collisions[i].particle2], lastEarliest);
                    collisions[i].time = newCollision.time;
                }
            } else {
                // the scanned collision is a particle-wall collision
                if (collisions[i].particle1 == particle1.id || collisions[i].particle1 == particle2.id) {
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
    for (int id = 0; id < N; id++) {
        if (id != particle1 && id != particle2) {
            particles[id].x = particles[id].x + (time - particles[id].tLast) * particles[id].vx;
            particles[id].y = particles[id].y + (time - particles[id].tLast) * particles[id].vy;
            particles[id].tLast = time;
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
    for (int id = 0; id < N; id++) {
        particles[id].x = particles[id].x + (endTime - particles[id].tLast) * particles[id].vx;
        particles[id].y = particles[id].y + (endTime - particles[id].tLast) * particles[id].vy;
        printf("%.6f, %.6f, %d, %d\n", particles[id].x, particles[id].y, particles[id].countW, particles[id].countP);
    }
    finalize();
    return 1;
}
