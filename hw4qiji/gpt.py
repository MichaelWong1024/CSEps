import sys

# Define the Particle and PotentialCollision classes
class Particle:
    def __init__(self, x, y, u, v, r):
        self.x = x
        self.y = y
        self.u = u
        self.v = v
        self.radius = r
        self.last_update = 0
        self.cwall = 0
        self.cpart = 0

class PotentialCollision:
    def __init__(self, particle1, particle2, collision_time, collision_type):
        self.particle1 = particle1
        self.particle2 = particle2
        self.collision_time = collision_time
        self.collision_type = collision_type # 'wall' or 'particle'

def update_position(particle, time):
    if time > particle.last_update:
        dt = time - particle.last_update
        particle.x += particle.u * dt
        particle.y += particle.v * dt
        particle.last_update = time

def calculate_wall_collision_time(particle, L_x, L_y):
    # Time for left/right wall
    time_x = ((L_x - particle.radius) if particle.u > 0 else particle.radius) - particle.x
    time_x /= abs(particle.u)
    # Time for top/bottom wall
    time_y = ((L_y - particle.radius) if particle.v > 0 else particle.radius) - particle.y
    time_y /= abs(particle.v)
    # Return the minimum positive time
    return min(t for t in [time_x, time_y] if t > 0)

def calculate_particle_collision_time(p1, p2):
    # Relative position and velocity
    dx = p2.x - p1.x
    dy = p2.y - p1.y
    dvx = p2.u - p1.u
    dvy = p2.v - p1.v

    # Quadratic coefficients: a*t^2 + b*t + c = 0
    a = dvx**2 + dvy**2
    b = 2 * (dx * dvx + dy * dvy)
    c = dx**2 + dy**2 - (4 * p1.radius**2)

    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        return float('inf')  # No collision if discriminant is negative

    # Quadratic formula to find the time t
    t1 = (-b + discriminant**0.5) / (2 * a)
    t2 = (-b - discriminant**0.5) / (2 * a)

    # We need the smallest positive time
    positive_times = [t for t in [t1, t2] if t > 0]
    return min(positive_times, default=float('inf'))

def process_collision(collision):
    if collision.collision_type == 'wall':
        particle = collision.particle1
        update_position(particle, collision.collision_time)
        # Wall collision: reverse the velocity component
        if particle.x <= particle.radius or particle.x >= L_x - particle.radius:
            particle.u = -particle.u
        if particle.y <= particle.radius or particle.y >= L_y - particle.radius:
            particle.v = -particle.v
        particle.cwall += 1
    else:
        # Particle collision: swap velocities
        p1, p2 = collision.particle1, collision.particle2
        update_position(p1, collision.collision_time)
        update_position(p2, collision.collision_time)
        p1.u, p2.u = p2.u, p1.u
        p1.v, p2.v = p2.v, p1.v
        p1.cpart += 1
        p2.cpart += 1

# Parse command line arguments
filename = sys.argv[1]
endtime = float(sys.argv[2])

# Read particle information from file
particles = []
with open(filename, 'r') as file:
    num_particles = int(file.readline().strip())
    radius = float(file.readline().strip())
    L_x, L_y = map(float, file.readline().strip().split())
    for _ in range(num_particles):
        x, y, u, v = map(float, file.readline().strip().split())
        particles.append(Particle(x, y, u, v, radius))

# Calculate all potential collision times
potential_collisions = []
for i, p1 in enumerate(particles):
    for j in range(i + 1, len(particles)):
        p2 = particles[j]
        collision_time = calculate_particle_collision_time(p1, p2)
        if collision_time != float('inf'):
            potential_collisions.append(PotentialCollision(p1, p2, collision_time, 'particle'))

    wall_collision_time = calculate_wall_collision_time(p1, L_x, L_y)
    potential_collisions.append(PotentialCollision(p1, None, wall_collision_time, 'wall'))

# Sort the potential collisions by time
potential_collisions.sort(key=lambda c: c.collision_time)

# Main simulation loop
current_time = 0
while potential_collisions and potential_collisions[0].collision_time <= endtime:
    # Sort the collisions array each time to ensure the earliest collision is processed next
    potential_collisions.sort(key=lambda c: c.collision_time)
    next_collision = potential_collisions.pop(0)
    
    if next_collision.collision_time > current_time:
        # Update current time to the time of collision
        current_time = next_collision.collision_time

        # Process the collision
        process_collision(next_collision)

        # Re-compute potential collisions for involved particles
        # If it was a wall collision, recompute for that one particle
        if next_collision.collision_type == 'wall':
            p = next_collision.particle1
            wall_collision_time = calculate_wall_collision_time(p, L_x, L_y)
            potential_collisions.append(PotentialCollision(p, None, current_time + wall_collision_time, 'wall'))
        else:  # It was a particle collision, recompute for both particles
            p1, p2 = next_collision.particle1, next_collision.particle2
            for p in [p1, p2]:
                wall_collision_time = calculate_wall_collision_time(p, L_x, L_y)
                potential_collisions.append(PotentialCollision(p, None, current_time + wall_collision_time, 'wall'))
        
        # Recompute particle-particle collision times
        for i, p1 in enumerate(particles):
            for j in range(i + 1, len(particles)):
                p2 = particles[j]
                # Only recompute for the particles involved in the last collision
                if p1 == next_collision.particle1 or p1 == next_collision.particle2 or \
                   p2 == next_collision.particle1 or p2 == next_collision.particle2:
                    collision_time = calculate_particle_collision_time(p1, p2)
                    if collision_time != float('inf'):
                        potential_collisions.append(PotentialCollision(p1, p2, current_time + collision_time, 'particle'))

# Update all particle positions to the end time
for particle in particles:
    update_position(particle, endtime)
    print(f"{particle.x:.6f}, {particle.y:.6f}, {particle.cwall}, {particle.cpart}")

# Exit the program
sys.exit(0)
