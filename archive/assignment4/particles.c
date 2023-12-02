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

/* Read simulation parameters and particle data from file */
void read(const char* fileName) {
    FILE* fp = fopen(fileName, "r"); // Open file for reading

    // Read simulation parameters
    fscanf(fp, "%d %lf %lf %lf", &numP, &radius, &lenX, &lenY);

    // Allocate memory for particles and collision arrays
    particles = malloc(numP * sizeof(Particle));
    collisions = malloc(numP * (numP + 1) / 2 * sizeof(Collision));

    // Populate particle array with data from file
    for (int i = 0; i < numP; i++) {
        fscanf(fp, "%lf %lf %lf %lf", 
               &particles[i].x, &particles[i].y, 
               &particles[i].vx, &particles[i].vy);
        particles[i].id = i;
        particles[i].tLast = 0;
        particles[i].countP = 0;
        particles[i].countW = 0;
    }

    fclose(fp); // Close file
}

/* Calculate the particle-wall collision based on the given particle */
Collision findColW(Particle p) {
    // Find the collision time with each wall
    double tLeft = (p.vx < 0) ? (p.tLast + (radius - p.x) / p.vx) : INFINITY;
    double tRight = (p.vx > 0) ? (p.tLast + ((lenX - radius) - p.x) / p.vx) : INFINITY;
    double tBottom = (p.vy < 0) ? (p.tLast + (radius - p.y) / p.vy) : INFINITY;
    double tTop = (p.vy > 0) ? (p.tLast + ((lenY - radius) - p.y) / p.vy) : INFINITY;

    // Determine the earliest collision time
    double t = minT(minT(tLeft, tRight), minT(tTop, tBottom));

    // Populate the collision structure
    Collision collision = {p.id, -1, t};  // -1 indicates a wall collision
    return collision;
}

/* Calculate the particle-particle collision based on the two given particles */
Collision findColP(Particle p1, Particle p2, Collision lastCol) {
    // Synchronize particle times, updating positions if needed
    double tDiff = p1.tLast - p2.tLast;
    if (fabs(tDiff) > 0.0000001) {
        Particle* pToUpdate = tDiff > 0 ? &p2 : &p1;
        tDiff = fabs(tDiff);
        pToUpdate->x += tDiff * pToUpdate->vx;
        pToUpdate->y += tDiff * pToUpdate->vy;
    }

    // Calculate collision time using quadratic formula components
    double dx = p1.x - p2.x, dvx = p1.vx - p2.vx;
    double dy = p1.y - p2.y, dvy = p1.vy - p2.vy;
    double a = dvx * dvx + dvy * dvy;
    double b = 2 * (dx * dvx + dy * dvy);
    double c = dx * dx + dy * dy - 4 * radius * radius;
    double discriminant = b * b - 4 * a * c;

    // Determine if and when the collision occurs
    double t = INFINITY;
    if (discriminant >= 0) {
        double sqrtDisc = sqrt(discriminant);
        double t1 = (-b - sqrtDisc) / (2 * a);
        t = t1 >= 0 ? t1 : INFINITY;  // Only consider future collisions
    }

    // Construct the collision result
    Collision collision = {p1.id, p2.id, t + lastCol.t};
    return collision;
}

/* Initialize the collision array */
void initColArr() {
    // Set the initial state for the last collision
    lastCol = (Collision){-1, -1, 0};

    // Initialize particle-wall collisions in the array
    for (int i = 0; i < numP; ++i) {
        collisions[i] = findColW(particles[i]);
    }
    
    // Initialize particle-particle collisions
    for (int i = 0, id = numP; i < numP; ++i) {
        for (int j = i + 1; j < numP; ++j) {
            collisions[id++] = findColP(particles[i], particles[j], lastCol);
        }
    }
}

/* Sort the collision array using insertion sort */
void sortColArr() {
    int totalCollisions = numP * (numP + 1) / 2; // Calculate the total number of collisions once

    // Iterate over each collision, starting from the second element
    for (int i = 1; i < totalCollisions; ++i) {
        Collision current = collisions[i];
        int j = i - 1;

        // Shift collisions that have greater time than current
        while (j >= 0 && collisions[j].t > current.t) {
            collisions[j + 1] = collisions[j];
            --j;
        }

        // Place current collision in its correct position
        collisions[j + 1] = current;
    }
}

/* Update the state of a particle after a wall collision */
void updatePW(Particle particle, Collision earliest) {
    // Calculate the new location of the particle
    Particle *p = &particles[particle.id]; // Direct reference to the particle in the array
    p->x += (earliest.t - p->tLast) * p->vx;
    p->y += (earliest.t - p->tLast) * p->vy;
    
    // Determine collision with vertical or horizontal walls and update velocity
    if (isEq(p->x, radius, 0.0000001f) || isEq(p->x, lenX - radius, 0.0000001f)) {
        p->vx = -p->vx; // Reflect horizontal velocity for vertical wall collision
    }
    if (isEq(p->y, radius, 0.0000001f) || isEq(p->y, lenY - radius, 0.0000001f)) {
        p->vy = -p->vy; // Reflect vertical velocity for horizontal wall collision
    }
    
    // Update the collision time and increment the wall collision counter
    p->tLast = earliest.t;
    p->countW++;
}

/* Update the states of two particles after a particle-particle collision */
void updatePP(Particle p1, Particle p2, Collision earliest) {
    // Direct references to the particles in the array
    Particle *part1 = &particles[p1.id];
    Particle *part2 = &particles[p2.id];
    
    // Update the locations of both particles
    part1->x += (earliest.t - part1->tLast) * part1->vx;
    part1->y += (earliest.t - part1->tLast) * part1->vy;
    part2->x += (earliest.t - part2->tLast) * part2->vx;
    part2->y += (earliest.t - part2->tLast) * part2->vy;
    
    // Swap velocities assuming a perfectly elastic collision
    double tempVx = part1->vx;
    double tempVy = part1->vy;
    part1->vx = part2->vx;
    part1->vy = part2->vy;
    part2->vx = tempVx;
    part2->vy = tempVy;
    
    // Update the collision time for both particles
    part1->tLast = part2->tLast = earliest.t;
    
    // Increment the particle-particle collision counters
    part1->countP++;
    part2->countP++;
}

/* Advances the simulation to the next collision event */
void jumpCol() {
    // Retrieve the earliest collision
    Collision earliest = collisions[0];
    lastCol = earliest; // Update the last collision record

    if (earliest.p2 == -1) {
        // If p2 is -1, it's a particle-wall collision
        updatePW(particles[earliest.p1], earliest);
    } else {
        // Otherwise, it's a particle-particle collision
        updatePP(particles[earliest.p1], particles[earliest.p2], earliest);
    }
}

/* Updates the collision array with new collisions for involved particles */
void updateColArr() {
    Collision earliest = lastCol;
    int involvedId1 = earliest.p1;
    int involvedId2 = earliest.p2;
    
    // Loop over the collision array to update collisions involving the particles from the latest collision
    for (int i = 0; i < numP * (numP + 1) / 2; i++) {
        // Check if the current collision involves one of the particles from the latest collision
        if (collisions[i].p1 == involvedId1 || collisions[i].p2 == involvedId1 || 
            collisions[i].p1 == involvedId2 || collisions[i].p2 == involvedId2) {
            
            // Update particle-wall collision if p2 is -1 or update particle-particle collision otherwise
            collisions[i] = (collisions[i].p2 == -1) ? 
                            findColW(particles[collisions[i].p1]) : 
                            findColP(particles[collisions[i].p1], particles[collisions[i].p2], lastCol);
        }
    }
}

/* Updates positions of all particles not involved in the last collision to the time of the last collision */
void updateParArr() {
    // Use the timestamp of the last collision
    double t = lastCol.t;
    
    // Iterate through all particles
    for (int id = 0; id < numP; id++) {
        // Skip particles involved in the last collision
        if (id != lastCol.p1 && id != lastCol.p2) {
            // Update the position of the particle based on its velocity and the time elapsed since its last update
            particles[id].x += (t - particles[id].tLast) * particles[id].vx;
            particles[id].y += (t - particles[id].tLast) * particles[id].vy;
            // Record the time of the last update to the particle
            particles[id].tLast = t;
        }
    }
}

/* free the allocated memory */
void finalize() {
    free(particles);
    free(collisions);
}
