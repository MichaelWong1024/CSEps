#include "particles.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Utility function to calculate the time to wall collision
double time_to_wall_collision(Particle *p, double box_x, double box_y, double radius) {
    double time_x, time_y;
    // Check for collision with left or right wall
    if (p->vx > 0) {
        time_x = (box_x - radius - p->x) / p->vx;
    } else {
        time_x = (radius - p->x) / p->vx;
    }
    // Check for collision with top or bottom wall
    if (p->vy > 0) {
        time_y = (box_y - radius - p->y) / p->vy;
    } else {
        time_y = (radius - p->y) / p->vy;
    }
    // Return the time for the earliest wall collision
    return fmin(time_x, time_y);
}

// Utility function to calculate the time to particle collision
double time_to_particle_collision(Particle *p1, Particle *p2, double radius) {
    // Vector from p1 to p2
    double dx = p2->x - p1->x;
    double dy = p2->y - p1->y;
    // Relative velocity vector
    double dvx = p2->vx - p1->vx;
    double dvy = p2->vy - p1->vy;

    // Coefficients of the quadratic equation (At^2 + Bt + C = 0)
    double A = dvx * dvx + dvy * dvy;
    double B = 2 * (dx * dvx + dy * dvy);
    double C = dx * dx + dy * dy - 4 * radius * radius;

    // Discriminant
    double D = B * B - 4 * A * C;
    if (D < 0) return DBL_MAX; // No collision if discriminant is negative

    double t1 = (-B - sqrt(D)) / (2 * A);
    double t2 = (-B + sqrt(D)) / (2 * A);

    // We only consider future collisions
    double t = (t1 > 0) ? t1 : ((t2 > 0) ? t2 : DBL_MAX);

    return t;
}

// Function to initialize particles from a file
void initialize_particles(Particle *particles, int num_particles, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Skip the first two lines which contain the number of particles and radius
    fscanf(file, "%*[^\n]\n");
    fscanf(file, "%*[^\n]\n");

    for (int i = 0; i < num_particles; i++) {
        fscanf(file, "%lf %lf %lf %lf\n", &particles[i].x, &particles[i].y, &particles[i].vx, &particles[i].vy);
        particles[i].last_updated_time = 0;
        particles[i].cwall = 0;
        particles[i].cpart = 0;
    }

    fclose(file);
}

// Function to calculate all potential collision times
void calculate_all_potential_collisions(Particle *particles, int num_particles, Collision *collisions, double box_x, double box_y, double radius) {
    int collision_index = 0;

    for (int i = 0; i < num_particles; i++) {
        // Calculate wall collisions for each particle
        collisions[collision_index].time = time_to_wall_collision(&particles[i], box_x, box_y, radius);
        collisions[collision_index].p1 = &particles[i];
        collisions[collision_index].p2 = NULL;
        collision_index++;

        // Calculate particle-particle collisions
        for (int j = i + 1; j < num_particles; j++) {
            collisions[collision_index].time = time_to_particle_collision(&particles[i], &particles[j], radius);
            collisions[collision_index].p1 = &particles[i];
            collisions[collision_index].p2 = &particles[j];
            collision_index++;
        }
    }
}

// Function to sort collisions by time
void sort_collisions(Collision *collisions, int num_collisions) {
    // Insertion sort algorithm
    for (int i = 1; i < num_collisions; i++) {
        Collision key = collisions[i];
        int j = i - 1;

        while (j >= 0 && collisions[j].time > key.time) {
            collisions[j + 1] = collisions[j];
            j--;
        }
        collisions[j + 1] = key;
    }
}

// Function to process a collision
void process_collision(Collision *collision, double box_x, double box_y, double radius) {
    Particle *p1 = collision->p1;
    Particle *p2 = collision->p2;
    double time = collision->time;

    // Update positions to the collision time
    update_position(p1, time);
    if (p2) {
        update_position(p2, time);
        particle_particle_collision(p1, p2, radius);
    } else {
        particle_wall_collision(p1, collision->wall, box_x, box_y, radius);
    }
}

// Function to update particle positions
void update_particle_positions(Particle *particles, int num_particles, double end_time) {
    for (int i = 0; i < num_particles; i++) {
        update_position(&particles[i], end_time);
    }
}

// Function to free memory
void free_memory(Particle *particles, Collision *collisions) {
    free(particles);
    free(collisions);
}
