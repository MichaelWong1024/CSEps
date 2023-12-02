import math
import itertools
import sys

# Class to represent a particle
class Particle:
    def __init__(self, x, y, vx, vy, radius):
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
        self.radius = radius
        self.last_update_time = 0.0
        self.wall_collision_count = 0
        self.particle_collision_count = 0

    def update_position(self, time):
        elapsed_time = time - self.last_update_time
        self.x += self.vx * elapsed_time
        self.y += self.vy * elapsed_time
        self.last_update_time = time

    def collide_with_wall(self, wall):
        if wall == 'vertical':
            self.vx = -self.vx
        else:
            self.vy = -self.vy
        self.wall_collision_count += 1

    def collide_with_particle(self, other):
        self.vx, other.vx = other.vx, self.vx
        self.vy, other.vy = other.vy, self.vy
        self.particle_collision_count += 1
        other.particle_collision_count += 1

# Class to represent a collision
class Collision:
    def __init__(self, time, particles, collision_type):
        self.time = time
        self.particles = particles
        self.collision_type = collision_type

# Class to simulate the particle dynamics
class Simulation:
    def __init__(self, particles, box_length_x, box_length_y):
        self.particles = particles
        self.box_length_x = box_length_x
        self.box_length_y = box_length_y
        self.collisions = []

    def calculate_wall_collisions(self, particle):
        times_to_collision = {'left': float('inf'), 'right': float('inf'),
                              'bottom': float('inf'), 'top': float('inf')}
        if particle.vx > 0:
            times_to_collision['right'] = (self.box_length_x - particle.radius - particle.x) / particle.vx
        elif particle.vx < 0:
            times_to_collision['left'] = (particle.radius - particle.x) / particle.vx
        if particle.vy > 0:
            times_to_collision['top'] = (self.box_length_y - particle.radius - particle.y) / particle.vy
        elif particle.vy < 0:
            times_to_collision['bottom'] = (particle.radius - particle.y) / particle.vy

        earliest_collision_time = min(times_to_collision.values())
        collision_wall = [wall for wall, time in times_to_collision.items() if time == earliest_collision_time][0]
        return Collision(earliest_collision_time + particle.last_update_time, [particle], collision_wall)

    def calculate_particle_collisions(self, particle, other):
        dx = other.x - particle.x
        dy = other.y - particle.y
        dvx = other.vx - particle.vx
        dvy = other.vy - particle.vy

        a = dvx**2 + dvy**2
        b = 2 * (dx * dvx + dy * dvy)
        c = dx**2 + dy**2 - (2 * particle.radius)**2

        discriminant = b**2 - 4 * a * c
        if discriminant < 0:
            return Collision(float('inf'), [particle, other], 'particle')

        discriminant_sqrt = math.sqrt(discriminant)
        t1 = (-b - discriminant_sqrt) / (2 * a)
        t2 = (-b + discriminant_sqrt) / (2 * a)

        if t1 > 0:
            collision_time = t1
        elif t2 > 0:
            collision_time = t2
        else:  # both times are in the past
            return Collision(float('inf'), [particle, other], 'particle')

        return Collision(collision_time + particle.last_update_time, [particle, other], 'particle')

    def simulate(self, end_time):
        # Initialize potential collisions
        for particle in self.particles:
            self.collisions.append(self.calculate_wall_collisions(particle))
        for particle, other in itertools.combinations(self.particles, 2):
            self.collisions.append(self.calculate_particle_collisions(particle, other))

        # Process collisions until end time
        while True:
            self.collisions.sort(key=lambda collision: collision.time)
            next_collision = self.collisions.pop(0)
            if next_collision.time > end_time:
                break
            for particle in next_collision.particles:
                particle.update_position(next_collision.time)
            if next_collision.collision_type == 'vertical' or next_collision.collision_type == 'horizontal':
                next_collision.particles[0].collide_with_wall(next_collision.collision_type)
            else:
                next_collision.particles[0].collide_with_particle(next_collision.particles[1])

            # Recalculate potential collisions involving the colliding particles
            self.collisions = [collision for collision in self.collisions if not any(p in collision.particles for p in next_collision.particles)]
            for particle in next_collision.particles:
                self.collisions.append(self.calculate_wall_collisions(particle))
                for other in [p for p in self.particles if p != particle]:
                    self.collisions.append(self.calculate_particle_collisions(particle, other))

        # Update all particle positions to the end time
        for particle in self.particles:
            particle.update_position(end_time)

    def get_simulation_results(self):
        return [(particle.x, particle.y, particle.wall_collision_count, particle.particle_collision_count) for particle in self.particles]

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 particle.py <filename> <end_time>")
        sys.exit(1)

    filename = sys.argv[1]
    end_time = float(sys.argv[2])
    
    # Read the file to get particle data
    with open(filename, 'r') as file:
        N = int(file.readline().strip())
        radius = float(file.readline().strip())
        box_length_x, box_length_y = map(float, file.readline().strip().split())
        particles = []
        for line in file:
            x, y, vx, vy = map(float, line.strip().split())
            particles.append(Particle(x, y, vx, vy, radius))

    # Run the simulation
    simulation = Simulation(particles, box_length_x, box_length_y)
    simulation.simulate(end_time)
    results = simulation.get_simulation_results()

    # Print the results
    for result in results:
        print("{:.6f}, {:.6f}, {}, {}".format(*result))
        
# Example usage:
# particles = [Particle(x, y, vx, vy, radius) for x, y, vx, vy in initial_conditions]
# simulation = Simulation(particles, box_length_x, box_length_y)
# simulation.simulate(end_time)
# results = simulation.get_simulation_results()
