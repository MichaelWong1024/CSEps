import sys
import math


class Particle:
    def __init__(self, xpos, ypos, xvel, yvel):
        self.xpos = xpos
        self.ypos = ypos
        self.xvel = xvel
        self.yvel = yvel
        self.last_update = 0.0
        self.cwall = 0
        self.cpart = 0

class Collision:
    def __init__(self, time, p1, p2, ctype):
        self.time = time
        self.p1 = p1
        self.p2 = p2
        self.ctype = ctype  # 'wall' or 'particle'

def read_particles_from_file(filename):
    particles = []
    with open(filename, 'r') as file:
        num_particles = int(file.readline().strip())
        radius = float(file.readline().strip())
        box_size = tuple(map(float, file.readline().strip().split()))
        
        for _ in range(num_particles):
            xpos, ypos, xvel, yvel = map(float, file.readline().strip().split())
            particles.append(Particle(xpos, ypos, xvel, yvel))
    return particles, radius, box_size


def calculate_wall_collision_times(particles, box_size, radius):
    collisions = []
    L_x, L_y = box_size
    for i, particle in enumerate(particles):
        # Calculate time to hit left/right walls
        if particle.xvel > 0:
            time_x = (L_x - radius - particle.xpos) / particle.xvel
            collisions.append(Collision(time_x, i, -1, 'wall'))
        elif particle.xvel < 0:
            time_x = (radius - particle.xpos) / particle.xvel
            collisions.append(Collision(time_x, i, -2, 'wall'))
        
        # Calculate time to hit top/bottom walls
        if particle.yvel > 0:
            time_y = (L_y - radius - particle.ypos) / particle.yvel
            collisions.append(Collision(time_y, i, -3, 'wall'))
        elif particle.yvel < 0:
            time_y = (radius - particle.ypos) / particle.yvel
            collisions.append(Collision(time_y, i, -4, 'wall'))
        
    return collisions


def calculate_particle_collision_times(particles, radius):
    collisions = []
    num_particles = len(particles)

    for i in range(num_particles):
        for j in range(i+1, num_particles):
            p1 = particles[i]
            p2 = particles[j]

            # Calculate the relative position and velocity vectors
            dx = p2.xpos - p1.xpos
            dy = p2.ypos - p1.ypos
            dvx = p2.xvel - p1.xvel
            dvy = p2.yvel - p1.yvel
            
            # Calculate the coefficients of the quadratic equation
            a = dvx**2 + dvy**2
            b = 2 * (dx * dvx + dy * dvy)
            c = dx**2 + dy**2 - (2*radius)**2
            
            # Calculate the discriminant
            discriminant = b**2 - 4*a*c
            
            # Only proceed if the discriminant is non-negative
            if discriminant >= 0:
                # Compute the two possible collision times
                t1 = (-b + math.sqrt(discriminant)) / (2*a)
                t2 = (-b - math.sqrt(discriminant)) / (2*a)
                # We are interested in the smallest positive time
                if t1 > 0 and t2 > 0:
                    collision_time = min(t1, t2)
                elif t1 > 0:
                    collision_time = t1
                elif t2 > 0:
                    collision_time = t2
                else:
                    # Both times are negative, no collision will occur
                    continue
                
                collisions.append(Collision(collision_time, i, j, 'particle'))
    
    return collisions


def sort_collisions(collisions):
    # Sort the collisions based on the time
    collisions.sort(key=lambda collision: collision.time)


def update_particle_positions(particles, time):
    for particle in particles:
        if time > particle.last_update:
            elapsed_time = time - particle.last_update
            particle.xpos += particle.xvel * elapsed_time
            particle.ypos += particle.yvel * elapsed_time
            particle.last_update = time


def process_collision(collision, particles, radius):
    # If it's a wall collision
    if collision.ctype == 'wall':
        particle = particles[collision.p1]
        if collision.p2 == -1:  # Right wall
            particle.xvel = -particle.xvel
            particle.xpos = L_x - radius  # Update position to the point of collision
        elif collision.p2 == -2:  # Left wall
            particle.xvel = -particle.xvel
            particle.xpos = radius  # Update position to the point of collision
        elif collision.p2 == -3:  # Top wall
            particle.yvel = -particle.yvel
            particle.ypos = L_y - radius  # Update position to the point of collision
        elif collision.p2 == -4:  # Bottom wall
            particle.yvel = -particle.yvel
            particle.ypos = radius  # Update position to the point of collision
        particle.cwall += 1  # Increment wall collision count

    # If it's a particle collision
    elif collision.ctype == 'particle':
        p1 = particles[collision.p1]
        p2 = particles[collision.p2]
        
        # Swap velocities
        p1.xvel, p2.xvel = p2.xvel, p1.xvel
        p1.yvel, p2.yvel = p2.yvel, p1.yvel

        # Update positions to the point of collision
        p1.xpos += p1.xvel * (collision.time - p1.last_update)
        p1.ypos += p1.yvel * (collision.time - p1.last_update)
        p2.xpos += p2.xvel * (collision.time - p2.last_update)
        p2.ypos += p2.yvel * (collision.time - p2.last_update)

        p1.last_update = collision.time
        p2.last_update = collision.time
        p1.cpart += 1  # Increment particle collision count for p1
        p2.cpart += 1  # Increment particle collision count for p2


def update_collisions(collisions, particles, box_size, radius, processed_collision):
    # Remove old collisions involving the particles that just collided
    collisions = [c for c in collisions if c.p1 not in (processed_collision.p1, processed_collision.p2) and c.p2 not in (processed_collision.p1, processed_collision.p2)]

    # Recalculate wall collisions for the particles that just collided
    if processed_collision.ctype == 'wall':
        collisions.extend(calculate_wall_collision_times([particles[processed_collision.p1]], box_size, radius))
    elif processed_collision.ctype == 'particle':
        collisions.extend(calculate_wall_collision_times([particles[processed_collision.p1], particles[processed_collision.p2]], box_size, radius))

    # Recalculate particle collisions for the particles that just collided
    involved_particles = [processed_collision.p1]
    if processed_collision.ctype == 'particle':
        involved_particles.append(processed_collision.p2)

    for i in involved_particles:
        for j in range(len(particles)):
            if j not in involved_particles:
                new_collision_time = calculate_particle_collision_time(particles[i], particles[j], radius)
                if new_collision_time is not None:
                    collisions.append(Collision(new_collision_time, i, j, 'particle'))

    # Sort collisions again after updating
    sort_collisions(collisions)
    return collisions


def print_final_positions(particles):
    for i, particle in enumerate(particles):
        print(f"particle {i}: ({particle.xpos:.6f}, {particle.ypos:.6f}) ({particle.xvel:.6f}, {particle.yvel:.6f}): {particle.last_update:.6f} : ({particle.cwall}, {particle.cpart})")


def main():
    if len(sys.argv) != 3:
        print("Usage: python particle.py <filename> <end_time>")
        return

    filename = sys.argv[1]
    end_time = float(sys.argv[2])
    box_size = (0, 0)  # This will be updated when reading the file

    particles = read_particles_from_file(filename)
    collisions = calculate_wall_collision_times(particles, box_size, radius) + calculate_particle_collision_times(particles, radius)

    current_time = 0.0
    while current_time < end_time:
        sort_collisions(collisions)
        next_collision = collisions[0]
        
        if next_collision.time > end_time:
            break

        update_particle_positions(particles, next_collision.time)
        process_collision(next_collision, particles)
        current_time = next_collision.time
        update_collisions(collisions, particles, box_size, next_collision)

    update_particle_positions(particles, end_time)
    print_final_positions(particles)

if __name__ == "__main__":
    main()
