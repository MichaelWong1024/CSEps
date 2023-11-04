#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
///@brief Structure of node in the graph.
typedef struct particle{
    int ID;
    double x_position;
    double y_position;
    double x_velocity;
    double y_velocity;
    double last_updated_time;
    int count_wall_collisions;
    int count_particle_collisions;
}Particle;

///@brief Structure of node in the graph.
typedef struct collision{
    Particle* particle1;
    Particle* particle2; // NULL for wall collisions
    double potential_collision_time;
    bool iswall; 
    bool isX;
}Collision;
///--------------------------------Particle-Particle Function-------------------------------------------------------------------------------
double x_path(Particle p, double t) {
    return p.x_position + p.x_velocity * (t - p.last_updated_time);
}

double y_path(Particle p, double t) {
    return p.y_position + p.y_velocity * (t - p.last_updated_time);
}

double particle_particle_potential_collision_time(Particle p1, Particle p2,double R,double current_time) {
    // Calculate the coefficients for the quadratic equation
    double p1_x_pos = x_path(p1,current_time); 
    double p2_x_pos = x_path(p2,current_time);
    double p1_y_pos = y_path(p1,current_time);  
    double p2_y_pos = y_path(p2,current_time); 
    double a = p1_x_pos - p2_x_pos;
    double b = p1.x_velocity - p2.x_velocity;
    double c = p1_y_pos - p2_y_pos;
    double d = p1.y_velocity - p2.y_velocity;
    

    double A = b * b + d * d;
    double B = 2 * (a * b + c * d);
    double C = a * a + c * c - 4 * R * R;

    double discriminant = B * B - 4 * A * C;
    if (discriminant < 0) return -1; // No collision

    // Using the quadratic formula to find the roots
    double t1 = (-B - sqrt(discriminant)) / (2 * A)+current_time;
    double t2 = (-B + sqrt(discriminant)) / (2 * A)+current_time;
    printf("t1: %lf, t2: %lf, current_time: %lf\n", t1,t2,current_time);
    return (t1 < t2) ? t1 : t2; // Return the earlier collision time
}

void update_velocities_after_collision(Particle* p1, Particle* p2) {
    double temp_u = p1->x_velocity;
    double temp_v = p1->y_velocity;

    p1->x_velocity = p2->x_velocity;
    p1->y_velocity = p2->y_velocity;

    p2->x_velocity = temp_u;
    p2->y_velocity = temp_v;
}
void particle_particle_handle_actual_collision(Collision potential_collision) {
    double t_c = potential_collision.potential_collision_time;
    Particle* p1 = potential_collision.particle1;
    Particle* p2 = potential_collision.particle2;
    p1->x_position = x_path(*p1, t_c);
    p1->y_position = y_path(*p1, t_c);
    p2->x_position = x_path(*p2, t_c);
    p2->y_position = y_path(*p2, t_c);
    update_velocities_after_collision(p1, p2);
    p1->last_updated_time = t_c;
    p2->last_updated_time = t_c;
    p1->count_particle_collisions+=1;
    p2->count_particle_collisions+=1;

}
///--------------------------------------------------------------------------------------------------------------------------------------

///--------------------------------Particle-Wall Function-------------------------------------------------------------------------------
double wall_collision_time_left(Particle p, double R, bool isX) {
    double velocity = isX ? p.x_velocity : p.y_velocity;
    double position = isX ? p.x_position : p.y_position;

    if (velocity >= 0) return -1; // No future collision with left wall
    return p.last_updated_time + ((R - position) / velocity);
}

double wall_collision_time_right(Particle p, double L, double R, bool isX) {
    double velocity = isX ? p.x_velocity : p.y_velocity;
    double position = isX ? p.x_position : p.y_position;

    if (velocity <= 0) return -1; // No future collision with right wall
    return p.last_updated_time + (L - R - position) / velocity;
}

double wall_potential_collision_time(Particle p, double L, double R, bool isX) {
    double t_left = wall_collision_time_left(p, R, isX);
    double t_right = wall_collision_time_right(p, L, R, isX);
    double velocity = isX ? p.x_velocity : p.y_velocity;

    if (velocity < 0) return t_left;
    if (velocity > 0) return t_right;
    return -1; // No potential collision
}


void wall_handle_real_collision(Collision potential_collision,double R,double L) {
    double t_potential = potential_collision.potential_collision_time;
    Particle* p = potential_collision.particle1;
    bool isX = potential_collision.isX;
    double* position = isX ? &(p->x_position) : &(p->y_position);
    double* velocity = isX ? &(p->x_velocity) : &(p->y_velocity);

    if (t_potential == wall_collision_time_left(*p, R, isX)) {
        *position = R;
        *velocity = -*velocity;
    } else if (t_potential == wall_collision_time_right(*p, L, R, isX)) {
        *position = L - R;
        *velocity = -*velocity;
    }
    position = isX ? &(p->y_position) : &(p->x_position);
    *position = isX ? (y_path(*p,t_potential)):(x_path(*p,t_potential));
    p->last_updated_time = t_potential;
    p->count_wall_collisions+=1;
}
///--------------------------------------------------------------------------------------------------------------------------------------

///--------------------------------Find Collision Function-------------------------------------------------------------------------------
int find_all_earliest_wall_collisions(Particle* particles, Collision* collisions, int totalNum, double L_x, double L_y, double R, int collision_index, double current_time, double end_time) {
    for (int i = 0; i < totalNum; i++) {
        Particle* current_particle = &particles[i];

        // Calculate potential collision time with left/right wall
        double t_left_right = wall_potential_collision_time(*current_particle, L_x, R, true);

        // Calculate potential collision time with bottom/top wall
        double t_bottom_top = wall_potential_collision_time(*current_particle, L_y, R, false);

        // Determine the earliest collision time
        double earliest_collision_time = DBL_MAX;
        bool isX = false;

        if (t_left_right != -1 && (t_bottom_top == -1 || t_left_right < t_bottom_top)) {
            earliest_collision_time = t_left_right;
            isX = true;
        } else if (t_bottom_top != -1) {
            earliest_collision_time = t_bottom_top;
        }
        Collision* current_collision = &collisions[collision_index];
        
        // Store the collision details
        current_collision->particle1 = current_particle;
        current_collision->particle2 = NULL; // NULL for wall collisions
        current_collision->potential_collision_time = (earliest_collision_time > current_time)? earliest_collision_time: DBL_MAX;
        current_collision->isX = isX;
        current_collision->iswall = true;
        
        collision_index++;
    }
    return collision_index;
}

int find_all_particle_particle_collisions(Particle* particles, Collision* collisions, int totalNum, double R, int collision_index, double current_time,double end_time) {
    // Loop through all pairs of particles
    for (int i = 0; i < totalNum - 1; i++) {
        for (int j = i + 1; j < totalNum; j++) {
            // Calculate the potential collision time between particles i and j
            double potential_time = particle_particle_potential_collision_time(particles[i], particles[j], R, current_time);
            double p1_x_pos = x_path(particles[i],potential_time); 
            double p2_x_pos = x_path(particles[j],potential_time);
            double p1_y_pos = y_path(particles[i],potential_time);  
            double p2_y_pos = y_path(particles[j],potential_time); 
            double dx = p2_x_pos - p1_x_pos;
            double dy = p2_y_pos - p1_y_pos;
            if(fabs(sqrt(dx * dx + dy * dy) - 2 * R) < 1e-10){
                printf("ID1: %d, ID2: %d, time: %lf\n",particles[i].ID,particles[j].ID,potential_time);
            }
            // printf("time: %lf", potential_time);
            // If there's a potential collision
            // Store the collision details in the collisions array
            collisions[collision_index].particle1 = &particles[i];
            collisions[collision_index].particle2 = &particles[j];
            collisions[collision_index].potential_collision_time = (potential_time>current_time)? potential_time:DBL_MAX;
            collisions[collision_index].iswall = false; // This is a particle-particle collision, not a wall collision
            collisions[collision_index].isX = false;

            // Increment the collision_index
            collision_index++;
        }
    }
    return collision_index;
}
///--------------------------------------------------------------------------------------------------------------------------------------

///-------------------------------Process Collision Related Function-------------------------------------------------------------------------------------------------------
void sort_potential_collisions(Collision* collisions, int collision_count) {
    for (int i = 1; i < collision_count; i++) {
        Collision key = collisions[i];
        int j = i - 1;

        // Move elements of collisions[0..i-1] that are greater than the key's potential_collision_time
        // to one position ahead of their current position
        while (j >= 0 && collisions[j].potential_collision_time > key.potential_collision_time) {
            collisions[j + 1] = collisions[j];
            j = j - 1;
        }
        collisions[j + 1] = key;
    }
}

void process_earliest_collision(Collision* collisions, double R, double L_x, double L_y) {
    if (collisions[0].iswall) {
        double L = collisions[0].isX ? L_x : L_y;
        wall_handle_real_collision(collisions[0], R, L);
    } else {
        // It's a particle-particle collision
        particle_particle_handle_actual_collision(collisions[0]);
    }

}

void update_potential_collision_times(Collision* collisions, int collision_count, double R, double L_x, double L_y, double current_time, double end_time) {
    if (!collisions || collision_count <= 0) return;  // Check if collisions is not NULL and has elements

    Collision* earliest_collision = &collisions[0];  // Assuming the collisions array is already sorted
    Particle* p1 = earliest_collision->particle1;
    Particle* p2 = earliest_collision->particle2;

    for (int i = 0; i < collision_count; i++) {
        Collision* current_collision = &collisions[i];
        // printf("updating: %d\n",i);
        if (current_collision->particle1 && current_collision->particle2) { // part_part
            if ((current_collision->particle1->ID == p1->ID || (p2 && current_collision->particle1->ID == p2->ID)) ||
                (current_collision->particle2->ID == p1->ID || (p2 && current_collision->particle2->ID == p2->ID))) {
                double t = particle_particle_potential_collision_time(*current_collision->particle1, *current_collision->particle2, R, current_time);
                double p1_x_pos = x_path(*current_collision->particle1,t); 
                double p2_x_pos = x_path(*current_collision->particle2,t);
                double p1_y_pos = y_path(*current_collision->particle1,t);  
                double p2_y_pos = y_path(*current_collision->particle2,t); 
                double dx = p2_x_pos - p1_x_pos;
                double dy = p2_y_pos - p1_y_pos;
                if(fabs(sqrt(dx * dx + dy * dy) - 2 * R) < 1e-10){
                    printf("ID1: %d, ID2: %d, time: %lf\n",current_collision->particle1->ID,current_collision->particle2->ID,t);
                }
                current_collision->potential_collision_time = (t > current_time) ? t : DBL_MAX;
            }
        } else { // wall
            if (current_collision->particle1 && ((p1 && current_collision->particle1->ID == p1->ID) || (p2 && current_collision->particle1->ID == p2->ID))){
                double t_left_right = wall_potential_collision_time(*current_collision->particle1, L_x, R, true);

                // Calculate potential collision time with bottom/top wall
                double t_bottom_top = wall_potential_collision_time(*current_collision->particle1, L_y, R, false);

                // Determine the earliest collision time
                double earliest_collision_time = DBL_MAX;
                current_collision->isX = false;

                if (t_left_right != -1 && (t_bottom_top == -1 || t_left_right < t_bottom_top)) {
                    earliest_collision_time = t_left_right;
                     current_collision->isX = true;
                } else if (t_bottom_top != -1) {
                    earliest_collision_time = t_bottom_top;
                }
                current_collision->potential_collision_time = (earliest_collision_time > current_time) ? earliest_collision_time : DBL_MAX;
            }
        }
        
    }
}



void output_final_positions(Particle* particles, int totalNum, double lastCollisionTime) {
    for (int i = 0; i < totalNum; i++) {
        // Advance particle positions to the time of the last collision
        particles[i].x_position = x_path(particles[i], lastCollisionTime);
        particles[i].y_position = y_path(particles[i], lastCollisionTime);
        
        // Print the final position and collision counts for each particle
        printf("%.6f, %.6f, %d, %d\n", 
               particles[i].x_position, 
               particles[i].y_position, 
               particles[i].count_wall_collisions, 
               particles[i].count_particle_collisions);
    }
}
///-------------------------------Print Function-------------------------------------------------------------------------------------------------------
void print_collision(Collision* collisions, int collision_index) {
    printf("Collisions Details:\n");
    printf("-------------------------------------------------\n");
    for (int i = 0; i < collision_index; i++) {
        printf("Collision %d:\n", i + 1);
        
        printf("Particle 1 ID: %d\n", collisions[i].particle1->ID);
        if(collisions[i].particle2 != NULL) {
            printf("Particle 2 ID: %d\n", collisions[i].particle2->ID);
        } else {
            printf("Particle 2: Wall\n");
        }
        
        printf("Potential Collision Time: %lf\n", collisions[i].potential_collision_time);
        
        if(collisions[i].iswall) {
            printf("Collision Type: Wall\n");
            printf("isX: %d\n", collisions[i].isX);
        } else {
            printf("Collision Type: Particle\n");
            printf("isX: %d\n", collisions[i].isX);

        }
        
        printf("-------------------------------------------------\n");
    }
}
void print_particles(Particle* particles, int totalNum) {
    printf("Particle Details:\n");
    printf("-------------------------------------------------\n");
    for (int i = 0; i < totalNum; i++) {
        printf("Particle %d:\n", i);
        printf("X Position: %lf\n", particles[i].x_position);
        printf("Y Position: %lf\n", particles[i].y_position);
        printf("X Velocity: %lf\n", particles[i].x_velocity);
        printf("Y Velocity: %lf\n", particles[i].y_velocity);
        printf("Last Updated Time: %lf\n", particles[i].last_updated_time);
        printf("Wall Collisions: %d\n", particles[i].count_wall_collisions);
        printf("Particle Collisions: %d\n", particles[i].count_particle_collisions);
        printf("-------------------------------------------------\n");
    }
}
///--------------------------------------------------------------------------------------------------------------------------------------

///-------------------------------Initializatoin Function-------------------------------------------------------------------------------------------------------
Particle createParticles(double temp_x_position, double temp_y_position, double temp_x_velocity, double temp_y_velocity,int ID){
    Particle new_particle;
    new_particle.ID = ID;
    new_particle.x_position = temp_x_position;
    new_particle.y_position = temp_y_position;
    new_particle.x_velocity = temp_x_velocity;
    new_particle.y_velocity = temp_y_velocity;
    new_particle.last_updated_time = 0.0;
    new_particle.count_wall_collisions = 0;
    new_particle.count_particle_collisions = 0;
    return new_particle;
}
Particle* readFile(FILE* file,int* totalNum,double* R,double* L_x,double* L_y){
    double temp_x_position;double temp_y_position;double temp_x_velocity;double temp_y_velocity;
    int collision_index = 0;
    fscanf(file, "%d", totalNum);
    fscanf(file, "%lf", R);
    fscanf(file, "%lf %lf", L_x, L_y);
    Particle* particles = (Particle*)malloc(*totalNum * sizeof(Particle));
    

    for(int i = 0;i<*totalNum;i++){
        fscanf(file, "%lf %lf %lf %lf", &temp_x_position, &temp_y_position, &temp_x_velocity, &temp_y_velocity);
        particles[i] = createParticles(temp_x_position,temp_y_position,temp_x_velocity,temp_y_velocity,i);
    }
    return particles;
}
Collision* make_collision_array(Collision* potential_collisions,int* totalNum){
    potential_collisions = (Collision*)malloc((*totalNum*(*totalNum+1)/2) * sizeof(Particle));
    return potential_collisions;
}
///--------------------------------------------------------------------------------------------------------------------------------------




/// @brief free all the allocated memory insdie the graph: List, Node etc.
/// @param the input graph 
void finalize(Particle* particles, Collision* potential_collisions){
    free(particles);
    free(potential_collisions);
}
