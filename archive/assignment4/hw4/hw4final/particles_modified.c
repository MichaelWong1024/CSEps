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

/* Calculate the collision of a particle with the walls */
Collision findWC(const Particle* p) {
    double tX = (p->vx > 0) ? (lenX - radius - p->x) / p->vx : (radius - p->x) / p->vx;
    double tY = (p->vy > 0) ? (lenY - radius - p->y) / p->vy : (radius - p->y) / p->vy;
    
    Collision col;
    col.p1 = p->id;
    col.p2 = -1;  // Indicating a wall collision
    col.t = p->tLast + minT((p->vx != 0) ? tX : INFINITY, (p->vy != 0) ? tY : INFINITY);

    return col;
}

/* Calculate the particle-particle collision based on the two given particles */
Collision findPC(Particle* p1, Particle* p2) {
    // Equalize the time of both particles if they are not the same
    double timeDiff = p1->tLast - p2->tLast;
    if (!isEq(timeDiff, 0.0, 1e-7)) {
        if (timeDiff > 0) {
            p2->x += timeDiff * p2->vx;
            p2->y += timeDiff * p2->vy;
        } else {
            p1->x -= timeDiff * p1->vx; // Note that timeDiff is negative here
            p1->y -= timeDiff * p1->vy;
        }
    }
    
    // Calculate discriminant for collision time
    double dx = p1->x - p2->x;
    double dvx = p1->vx - p2->vx;
    double dy = p1->y - p2->y;
    double dvy = p1->vy - p2->vy;
    double discriminant = pow(2 * (dx * dvx + dy * dvy), 2) - 4 * (dvx * dvx + dvy * dvy) * (dx * dx + dy * dy - 4 * radius * radius);
    
    double time = INFINITY;
    if (discriminant >= 0) {
        double root = (-2 * (dx * dvx + dy * dvy) - sqrt(discriminant)) / (2 * (dvx * dvx + dvy * dvy));
        if (root >= 0) {
            time = root + minT(p1->tLast, p2->tLast);
        }
    }
    
    return (Collision){ .p1 = p1->id, .p2 = p2->id, .t = time };
}

/* Compute all initial collisions */
void computeAll(void) {
    // Initialize the lastEarliest global variable
    lastEarliest = (Collision){ .p1 = -1, .p2 = -1, .t = 0 };
    
    // Calculate particle-wall collisions for each particle
    for (int i = 0; i < nParticles; i++) {
        collisions[i] = findWC(&particles[i]);
    }
    
    // Calculate particle-particle collisions
    int index = nParticles; // Start index for particle-particle collisions
    for (int i = 0; i < nParticles; i++) {
        for (int j = i + 1; j < nParticles; j++) {
            collisions[index++] = findPC(&particles[i], &particles[j]);
        }
    }
}

/* Sort collisions using insertion sort */
void sortCols(void) {
    int totalCollisions = nParticles * (nParticles + 1) / 2;
    for (int i = 1; i < totalCollisions; i++) {
        Collision temp = collisions[i];
        int j = i - 1;
        while (j >= 0 && collisions[j].t > temp.t) {
            collisions[j + 1] = collisions[j];
            j--;
        }
        collisions[j + 1] = temp;
    }
}

/* Update particle after wall collision occurs */
void updateWC(Particle* p, const Collision* col) {
    // Update the particle's position based on the collision time
    p->x += (col->t - p->tLast) * p->vx;
    p->y += (col->t - p->tLast) * p->vy;

    // Reflect the velocity component if collision is with top/bottom walls
    if (isEq(p->x, radius, 0.0000001) || isEq(p->x, lenX - radius, 0.0000001)) {
        p->vx = -p->vx;
    }

    // Reflect the velocity component if collision is with left/right walls
    if (isEq(p->y, radius, 0.0000001) || isEq(p->y, lenY - radius, 0.0000001)) {
        p->vy = -p->vy;
    }

    // Update the last collision time to the time of this collision
    p->tLast = col->t;

    // Increment the wall collision count for the particle
    p->colWall++;
}

/* Update particles after particle-particle collision occurs */
void updatePC(Particle* p1, Particle* p2, const Collision* col) {
    // Update the location of particle1 based on the time of collision
    p1->x += (col->t - p1->tLast) * p1->vx;
    p1->y += (col->t - p1->tLast) * p1->vy;
    
    // Update the location of particle2 based on the time of collision
    p2->x += (col->t - p2->tLast) * p2->vx;
    p2->y += (col->t - p2->tLast) * p2->vy;
    
    // Swap the velocities of the two particles
    double tempVx = p1->vx;
    double tempVy = p1->vy;
    p1->vx = p2->vx;
    p1->vy = p2->vy;
    p2->vx = tempVx;
    p2->vy = tempVy;
    
    // Update the last collision time for both particles
    p1->tLast = p2->tLast = col->t;
    
    // Increment the particle-particle collision count for both particles
    p1->colPart++;
    p2->colPart++;
}

/* Process the next earliest collision */
void processNext(void) {
    Collision earliest = collisions[0];
    lastEarliest = earliest;  // Update the last earliest collision to the current one.

    if (earliest.p2 == -1) {
        // If p2 is -1, it's a particle-wall collision.
        updateWC(&particles[earliest.p1], &earliest);
    } else {
        // Otherwise, it's a particle-particle collision.
        updatePC(&particles[earliest.p1], &particles[earliest.p2], &earliest);
    }

    updateCols();  // Assuming updateCols will refresh the collisions array.
}

void updateCols() {
    Collision earliest = lastEarliest;
    int p1 = earliest.p1;
    int p2 = earliest.p2;

    // Iterate through all possible collision pairs
    for (int i = 0; i < nParticles * (nParticles - 1) / 2; i++) {
        // Check if the current collision involves one of the particles from the earliest collision
        if (collisions[i].p1 == p1 || collisions[i].p1 == p2 || 
            collisions[i].p2 == p1 || collisions[i].p2 == p2) {
            if (collisions[i].p2 == -1) { // Particle-wall collision
                collisions[i].t = findWC(&particles[collisions[i].p1]).t;
            } else { // Particle-particle collision
                int otherIndex = (collisions[i].p1 == p1 || collisions[i].p1 == p2) ? 
                                  collisions[i].p2 : collisions[i].p1;
                collisions[i].t = findPC(&particles[collisions[i].p1], &particles[otherIndex]).t;
            }
        }
    }
}

void updateAll() {
    double time = lastEarliest.t; // Time of the last collision
    for (int i = 0; i < nParticles; i++) {
        double deltaTime = time - particles[i].tLast; // Time elapsed since the last update of the particle
        particles[i].x += deltaTime * particles[i].vx; // Update x position
        particles[i].y += deltaTime * particles[i].vy; // Update y position
        particles[i].tLast = time; // Set the last update time to the time of the last collision
    }
}

/* free the allocated memory */
void finish(void) {
    free(particles);
    free(collisions);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("Usage: %s <filename> <endtime>\n", argv[0]);
        return 1;
    }
    
    char* filename = argv[1];
    double endTime = atof(argv[2]);
    double deltaTime;

    // Initialize the particle array and the collision array
    read(filename);
    computeAll();
    sortCols();

    // Iterate until the simulation reaches the end time
    while (collisions[0].t <= endTime) {
        processNext();
        updateCols();
        sortCols();
    }

    updateAll();
    // Update the positions of all the particles to the end time
    for (int index = 0; index < nParticles; index++) {
        deltaTime = endTime - particles[index].tLast;
        particles[index].x += deltaTime * particles[index].vx;
        particles[index].y += deltaTime * particles[index].vy;
        particles[index].tLast = endTime;
        printf("%.6f, %.6f, %d, %d\n",
               particles[index].x,
               particles[index].y,
               particles[index].colWall,
               particles[index].colPart);
    }

    finish();
    return 0;
}
