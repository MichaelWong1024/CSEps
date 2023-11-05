#include <stdio.h>
#include <stdlib.h>
#include "particles.h"

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <filename> <end_time>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *filename = argv[1];
    double end_time = atof(argv[2]);

    // Read the initial particle information from the file
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return EXIT_FAILURE;
    }

    int num_particles;
    double radius, box_x, box_y;
    fscanf(file, "%d", &num_particles);
    fscanf(file, "%lf", &radius);
    fscanf(file, "%lf %lf", &box_x, &box_y);
    fclose(file);

    // Allocate memory for particles and collisions
    Particle *particles = malloc(num_particles * sizeof(Particle));
    int num_collisions = num_particles * (num_particles + 1) / 2;
    Collision *collisions = malloc(num_collisions * sizeof(Collision));

    if (!particles || !collisions) {
        perror("Error allocating memory");
        return EXIT_FAILURE;
    }

    // Initialize particles and calculate all potential collisions
    initialize_particles(particles, num_particles, filename);
    calculate_all_potential_collisions(particles, num_particles, collisions, box_x, box_y, radius);

    // Main simulation loop
    double current_time = 0.0;
    while (current_time < end_time) {
        // Sort collisions
        sort_collisions(collisions, num_collisions);

        // Process earliest collision
        process_collision(&collisions[0]);
        current_time = collisions[0].time;

        // Update potential collisions
        for (int i = 0; i < num_collisions; i++) {
            if (collisions[i].p1 == collisions[0].p1 || collisions[i].p1 == collisions[0].p2 ||
                (collisions[i].p2 && (collisions[i].p2 == collisions[0].p1 || collisions[i].p2 == collisions[0].p2))) {
                // Recalculate collision time
                if (collisions[i].p2 == NULL) { // Wall collision
                    collisions[i].time = time_to_wall_collision(collisions[i].p1, radius, 'x');
                } else { // Particle collision
                    collisions[i].time = time_to_particle_collision(collisions[i].p1, collisions[i].p2, radius);
                }
            }
        }

        // Check if the next collision is beyond the end time
        if (collisions[0].time > end_time) {
            break;
        }
    }

    // Update all particles to the end time
    update_particle_positions(particles, num_particles, end_time);

    // Output the final positions and collision counts
    for (int i = 0; i < num_particles; ++i) {
        printf("%.6f, %.6f, %d, %d\n", particles[i].x, particles[i].y, particles[i].cwall, particles[i].cpart);
    }

    // Free allocated memory
    free_memory(particles, collisions);

    return EXIT_SUCCESS;
}

