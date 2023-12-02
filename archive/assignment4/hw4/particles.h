#include <stdbool.h>
typedef struct particle {
    int ID;
    double x_position;
    double y_position;
    double x_velocity;
    double y_velocity;
    double last_updated_time;
    int count_wall_collisions;
    int count_particle_collisions;
} Particle;

typedef struct collision {
    Particle* particle1;
    Particle* particle2; // NULL for wall collisions
    double potential_collision_time;
    bool iswall; 
    bool isX;
} Collision;

double x_path(Particle p, double t);
double y_path(Particle p, double t);
double particle_particle_potential_collision_time(Particle p1, Particle p2, double R);
void update_velocities_after_collision(Particle* p1, Particle* p2);
void particle_particle_handle_actual_collision(Collision potential_collision);
double wall_collision_time_left(Particle p, double R, bool isX);
double wall_collision_time_right(Particle p, double L, double R, bool isX);
double wall_potential_collision_time(Particle p, double L, double R, bool isX);
void wall_handle_real_collision(Collision potential_collision, double R, double L);
Particle createParticles(double temp_x_position, double temp_y_position, double temp_x_velocity, double temp_y_velocity, int ID);
int find_all_earliest_wall_collisions(Particle* particles, Collision* collisions, int totalNum, double L_x, double L_y, double R, int collision_index,double current_time, double end_time);
int find_all_particle_particle_collisions(Particle* particles, Collision* collisions, int totalNum, double R, int collision_index,double current_time,double end_time);
void sort_potential_collisions(Collision* collisions, int collision_count);
void process_earliest_collision(Collision* collisions, double R, double L_x, double L_y);
void update_potential_collision_times(Collision* collisions, int collision_count, double R, double L_x, double L_y,double current_time,double end_time);
void output_final_positions(Particle* particles, int totalNum, double lastCollisionTime);
void print_collision(Collision* collisions, int collision_index);
void print_particles(Particle* particles, int totalNum);
Particle* readFile(FILE* file, int* totalNum, double* R, double* L_x, double* L_y);
Collision* make_collision_array(Collision* potential_collisions, int* totalNum);
void finalize(Particle* particles, Collision* potential_collisions);
