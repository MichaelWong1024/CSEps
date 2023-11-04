// particles.h
#ifndef PARTICLES_H
#define PARTICLES_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define INFINITY 1.0e10

// Structure to represent a particle
typedef struct {
    double x, y;           // Position
    double vx, vy;         // Velocity
    double radius;         // Radius of the particle
    int wall_collisions;   // Count of wall collisions
    int particle_collisions; // Count of particle collisions
} Particle;

// Structure to represent a potential collision
typedef struct {
    int particleA, particleB;  // Indices of the particles involved (-1 for wall)
    double collision_time;     // Predicted time of the collision
} CollisionEvent;

// Function declarations
void simulate(Particle* particles, CollisionEvent* events, int num_particles, double box_length_x, double box_length_y, double end_time);
void update_collision_events(CollisionEvent* events, Particle* particles, int num_particles, double box_length_x, double box_length_y, int updated_particleA, int updated_particleB, double current_time);
int compare_collision_events(const void* a, const void* b);
double time_to_wall_collision(Particle* p, double box_length_x, double box_length_y);
double time_to_particle_collision(Particle* p1, Particle* p2, double radius);
void handle_wall_collision(Particle* p, double box_length_x, double box_length_y);
void handle_particle_collision(Particle* p1, Particle* p2);

#endif // PARTICLES_H
