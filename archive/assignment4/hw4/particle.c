#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Define a structure to represent a particle
typedef struct {
    double x, y;      // Positions
    double u, v;      // Velocities
    double last_updated; // Last update time
    int cwall, cpart; // Count of collisions with walls and other particles
} Particle;

// Define a structure to represent a potential collision
typedef struct {
    int particle1, particle2; // Indices of the particles (-1 for walls)
    double time;              // Collision time
} Collision;

// Function to calculate the time until collision with a wall
double wall_collision_time(Particle *p, double L_x, double L_y, double R) {
    double time_x = (p->u > 0) ? (L_x - R - p->x) / p->u : (R - p->x) / p->u;
    double time_y = (p->v > 0) ? (L_y - R - p->y) / p->v : (R - p->y) / p->v;
    if (p->u < 0) time_x = -time_x;
    if (p->v < 0) time_y = -time_y;
    return (time_x < time_y) ? time_x : time_y; // Return the smaller time
}

// Function to calculate the time until collision between two particles
double particle_collision_time(Particle *p1, Particle *p2, double R) {
    double dx = p1->x - p2->x;
    double dy = p1->y - p2->y;
    double dvx = p1->u - p2->u;
    double dvy = p1->v - p2->v;
    
    double a = dvx*dvx + dvy*dvy;
    double b = 2 * (dx*dvx + dy*dvy);
    double c = dx*dx + dy*dy - 4*R*R;
    
    double disc = b*b - 4*a*c;
    if (disc < 0) return INFINITY; // No collision if discriminant is negative
    
    double t1 = (-b - sqrt(disc)) / (2*a);
    double t2 = (-b + sqrt(disc)) / (2*a);
    
    return (t1 < t2) ? t1 : t2; // Return the smaller time
}

// Function to update the position of a particle
void update_position(Particle *p, double time) {
    if (time > p->last_updated) {
        double dt = time - p->last_updated;
        p->x += p->u * dt;
        p->y += p->v * dt;
        p->last_updated = time;
    }
}

// Function to process a wall collision
void process_wall_collision(Particle *p, double L_x, double L_y, double R) {
    // Check which wall the collision is with and update velocity accordingly
    if (p->x <= R || p->x >= L_x - R) {
        p->u = -p->u;
        p->cwall++;
    }
    if (p->y <= R || p->y >= L_y - R) {
        p->v = -p->v;
        p->cwall++;
    }
}

// Function to process a particle collision
void process_particle_collision(Particle *p1, Particle *p2) {
    // Exchange velocities
    double temp_u = p1->u;
    double temp_v = p1->v;
    p1->u = p2->u;
    p1->v = p2->v;
    p2->u = temp_u;
    p2->v = temp_v;
    
    // Update collision counts
    p1->cpart++;
    p2->cpart++;
}

// Function to read particles from a file
Particle *read_particles(const char *filename, int *N, double *R, double *L_x, double *L_y) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file.\n");
        exit(EXIT_FAILURE);
    }
    
    fscanf(file, "%d", N); // Number of particles
    fscanf(file, "%lf", R); // Radius
    fscanf(file, "%lf %lf", L_x, L_y); // Box dimensions
    
    Particle *particles = malloc(*N * sizeof(Particle));
    if (!particles) {
        fprintf(stderr, "Error allocating memory.\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < *N; ++i) {
        Particle *p = &particles[i];
        fscanf(file, "%lf %lf %lf %lf", &p->x, &p->y, &p->u, &p->v);
        p->last_updated = 0;
        p->cwall = p->cpart = 0;
    }
    
    fclose(file);
    return particles;
}

// Function to sort collisions by time
void sort_collisions(Collision *collisions, int count) {
    // Simple insertion sort
    for (int i = 1; i < count; ++i) {
        Collision temp = collisions[i];
        int j = i - 1;
        while (j >= 0 && collisions[j].time > temp.time) {
            collisions[j + 1] = collisions[j];
            j--;
        }
        collisions[j + 1] = temp;
    }
}

// Main simulation function
void simulate(Particle *particles, int N, double L_x, double L_y, double R, double end_time) {
    // Calculate the number of potential collisions
    int collision_count = N + (N*(N-1))/2;
    Collision *collisions = malloc(collision_count * sizeof(Collision));
    if (!collisions) {
        fprintf(stderr, "Error allocating memory for collisions.\n");
        free(particles);
        exit(EXIT_FAILURE);
    }
    
    // Initial calculation of potential collision times
    for (int i = 0; i < N; ++i) {
        // Calculate wall collision time
        collisions[i].particle1 = i;
        collisions[i].particle2 = -1; // Indicator for wall collision
        collisions[i].time = particles[i].last_updated + wall_collision_time(&particles[i], L_x, L_y, R);
        
        // Calculate particle-particle collision times
        for (int j = i + 1; j < N; ++j) {
            int index = N + i * (N - 1) / 2 + j - (i + 1);
            collisions[index].particle1 = i;
            collisions[index].particle2 = j;
            collisions[index].time = particles[i].last_updated + particle_collision_time(&particles[i], &particles[j], R);
        }
    }
    
    // Sort potential collisions by time
    sort_collisions(collisions, collision_count);
    
    double current_time = 0.0;
    while (current_time < end_time && collision_count > 0) {
        // Process the earliest potential collision
        Collision collision = collisions[0];
        current_time = collision.time;
        if (current_time > end_time) break; // Stop if we've reached the end time
        
        // Update particle positions to the collision time
        for (int i = 0; i < N; ++i) {
            update_position(&particles[i], current_time);
        }
        
        if (collision.particle2 == -1) { // Wall collision
            process_wall_collision(&particles[collision.particle1], L_x, L_y, R);
        } else { // Particle collision
            process_particle_collision(&particles[collision.particle1], &particles[collision.particle2]);
        }
        
        // Re-calculate potential collisions involving the colliding particles
        for (int i = 0; i < N; ++i) {
            if (i == collision.particle1 || i == collision.particle2) {
                collisions[i].time = current_time + wall_collision_time(&particles[i], L_x, L_y, R);
                for (int j = 0; j < N; ++j) {
                    if (j != i) {
                        int index = N + i * (N - 1) / 2 + j - (i + 1);
                        collisions[index].time = current_time + particle_collision_time(&particles[i], &particles[j], R);
                    }
                }
            }
        }
        
        // Sort the array of potential collisions again
        sort_collisions(collisions, collision_count);
    }
    
    // Advance all particle positions to the end time and output their final states
    for (int i = 0; i < N; ++i) {
        update_position(&particles[i], end_time);
        printf("%.6f, %.6f, %d, %d\n", particles[i].x, particles[i].y, particles[i].cwall, particles[i].cpart);
    }
    
    // Free allocated memory
    free(collisions);
    free(particles);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <filename> <end_time>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *filename = argv[1];
    double end_time = atof(argv[2]);

    if (end_time <= 0) {
        fprintf(stderr, "End time must be positive.\n");
        return EXIT_FAILURE;
    }

    int N; // Number of particles
    double R, L_x, L_y; // Radius and box dimensions

    // Read particles from the file
    Particle *particles = read_particles(filename, &N, &R, &L_x, &L_y);

    // Perform the simulation
    simulate(particles, N, L_x, L_y, R, end_time);

    return EXIT_SUCCESS;
}