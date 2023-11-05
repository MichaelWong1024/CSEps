#ifndef PARTICLES_H
#define PARTICLES_H

#include <stdbool.h>

#define MAX_COLLISIONS 1000  // Adjust based on expected maximum number of collisions

// Define the structure for a particle
typedef struct {
    double x, y;        // Position
    double vx, vy;      // Velocity
    double last_updated_time;
    int cwall;          // Number of collisions with walls
    int cpart;          // Number of collisions with other particles
} Particle;

// Define the structure for a potential collision
typedef struct {
    double time;         // Collision time
    Particle *p1;        // Pointer to the first particle involved in the collision
    Particle *p2;        // Pointer to the second particle involved (NULL for wall collisions)
} Collision;

// Function prototypes
void initialize_particles(Particle *particles, int num_particles, const char *filename);
void calculate_all_potential_collisions(Particle *particles, int num_particles, Collision *collisions, double box_x, double box_y, double radius);
void sort_collisions(Collision *collisions, int num_collisions);
void process_collision(Collision *collision);
void update_particle_positions(Particle *particles, int num_particles, double current_time);
void free_memory(Particle *particles, Collision *collisions);

#endif // PARTICLES_H

