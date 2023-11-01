#include <stdio.h>
#include <math.h>

// Define a structure to represent a particle
struct Particle {
    double x, y, u, v;
    double lastUpdateTime;
};

// Function to calculate the time of potential collision between two particles
double calculateCollisionTime(struct Particle p1, struct Particle p2, double R, double currentTime) {
    double deltax = p2.x - p1.x;
    double deltay = p2.y - p1.y;
    double deltavx = p2.u - p1.u;
    double deltavy = p2.v - p1.v;

    double A = pow(deltavx, 2) + pow(deltavy, 2);
    double B = 2 * (deltax * deltavx + deltay * deltavy);
    double C = pow(deltax, 2) + pow(deltay, 2) - pow(2 * R, 2);

    double discriminant = pow(B, 2) - 4 * A * C;
    if (discriminant < 0) return -1; // No collision

    double t1 = (-B - sqrt(discriminant)) / (2 * A);
    double t2 = (-B + sqrt(discriminant)) / (2 * A);

    // Choosing the earlier time after the current time
    if (t1 > currentTime) return t1;
    if (t2 > currentTime) return t2;

    return -1; // No future collision
}

// Function to update the particles' positions and velocities after a collision
void updateParticlesAfterCollision(struct Particle* p1, struct Particle* p2, double collisionTime) {
    // Updating positions at the collision time
    p1->x += p1->u * (collisionTime - p1->lastUpdateTime);
    p1->y += p1->v * (collisionTime - p1->lastUpdateTime);
    p2->x += p2->u * (collisionTime - p2->lastUpdateTime);
    p2->y += p2->v * (collisionTime - p2->lastUpdateTime);

    // Switching velocities (perfectly elastic collision, equal mass)
    double temp_u = p1->u, temp_v = p1->v;
    p1->u = p2->u;
    p1->v = p2->v;
    p2->u = temp_u;
    p2->v = temp_v;

    p1->lastUpdateTime = p2->lastUpdateTime = collisionTime;
}

// Main simulation loop where you'll call the above functions
int main() {
    // Initialization of particles and other parameters goes here

    // Main simulation loop: calculate collision times, sort, process collisions, update particles

    return 0;
}

