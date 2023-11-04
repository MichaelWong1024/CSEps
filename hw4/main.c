#include <stdio.h>
#include <stdlib.h>
#include "particles.h"

int main(int argc, char *argv[]) {
    FILE *file; 
    int totalNum; int collision_index= 0; int processed = 0;
    double R, L_x, L_y; 
    Collision* potential_collisions;
    Particle* particles;

    if (argc != 3) {
        printf("Usage: %s <file_name> <end_Time>\n", argv[0]);
        return 1;
    }

    // Extract the filename and particles from command-line arguments
    char *filename = argv[1];
    double endTime = atof(argv[2]);double currentTime = 0.0;

    // read file
    file = fopen(filename, "r");
    //avoid crash
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }
    else {
        // initialize the graph
        particles = readFile(file, &totalNum, &R, &L_x, &L_y);
        potential_collisions = make_collision_array(potential_collisions,&totalNum);
        fclose(file);
    }

    collision_index = find_all_earliest_wall_collisions(particles, potential_collisions, totalNum, L_x, L_y, R,collision_index,currentTime,endTime);
    collision_index = find_all_particle_particle_collisions(particles, potential_collisions, totalNum, R,collision_index,currentTime,endTime);

    sort_potential_collisions(potential_collisions, collision_index);
    print_collision(potential_collisions,collision_index);
    printf("\n\n\n");

while (currentTime < endTime  && potential_collisions[0].potential_collision_time <= endTime) {
    currentTime = potential_collisions[0].potential_collision_time;    
    // Process the earliest potential collision and update the current time
    process_earliest_collision(potential_collisions, R, L_x, L_y);
    // Update potential collision times for particles involved in the most recent collision
    update_potential_collision_times(potential_collisions, collision_index, R, L_x, L_y, currentTime, endTime);
    // Sort potential collisions again after updating
    sort_potential_collisions(potential_collisions, collision_index);
    // print_particles(particles,totalNum);

    // print_collision(potential_collisions, collision_index);
    processed++;
    printf("currentTime: %lf\n\n\n", currentTime);
}

    printf("currentTime: %lf\n\n\n", currentTime);
    // printf("currentTime: %lf\n", currentTime);
    print_collision(potential_collisions,collision_index);
    printf("Processed: %d\n",processed);
    output_final_positions(particles, totalNum, currentTime);
    finalize(particles, potential_collisions);
    return 0;
}