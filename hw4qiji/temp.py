class Particle:
    def __init__(self, id, xpos, ypos, xvel, yvel):
        self.id = id
        self.xpos = xpos
        self.ypos = ypos
        self.xvel = xvel
        self.yvel = yvel
        self.last_update = 0.0  # Time of the last update
        self.cwall = 0  # Count of collisions with walls
        self.cpart = 0  # Count of collisions with other particles

# Initialize an empty list to hold particles
particles = []

# Read the file and parse the data
with open(file_path, 'r') as file:
    lines = file.readlines()

# Extract the particle information
for i, line in enumerate(lines[3:]):  # Skip the first three lines as they contain other information
    parts = line.strip().split()
    xpos, ypos, xvel, yvel = map(float, parts)
    particle = Particle(i, xpos, ypos, xvel, yvel)
    particles.append(particle)

# Verify the initial state of each particle
for p in particles:
    print(f"particle {p.id}: ({p.xpos:.6f}, {p.ypos:.6f}) ({p.xvel:.6f}, {p.yvel:.6f}): {p.last_update:.6f} : ({p.cwall}, {p.cpart})")


import sys
from itertools import combinations

# Define constants for the wall identifiers
LEFT_WALL, RIGHT_WALL, LOWER_WALL, UPPER_WALL = 0, 1, 2, 3

# Define a large value to represent an impossible collision (or MAX time)
MAX_TIME = sys.float_info.max

# Define the collision structure
class Collision:
    def __init__(self, time, part1, part2, ctype):
        self.time = time  # Time at which collision occurs
        self.part1 = part1  # First particle involved in collision or wall
        self.part2 = part2  # Second particle involved in collision (or None for wall collisions)
        self.ctype = ctype  # Type of collision ("wall" or "particle")

# Assume the box's right and upper edge positions are provided in the file (L_x, L_y)
# Read the first three lines to get the number of particles, particle radius, and box dimensions
with open(file_path, 'r') as file:
    n_particles = int(file.readline().strip())
    radius = float(file.readline().strip())
    L_x, L_y = map(float, file.readline().strip().split())

# Function to calculate wall collision time
def calculate_wall_collision(particle, L, radius, velocity, is_vertical):
    if velocity > 0:
        wall_pos = L - radius
        wall = UPPER_WALL if is_vertical else RIGHT_WALL
    else:
        wall_pos = radius
        wall = LOWER_WALL if is_vertical else LEFT_WALL

    # Time = Distance / Velocity, considering the direction of the velocity
    time = (wall_pos - (particle.ypos if is_vertical else particle.xpos)) / velocity
    return time, wall

# Function to initialize the collision array
def initialize_collision_array(particles, L_x, L_y, radius):
    collision_array = []

    # Calculate potential wall collisions for each particle
    for p in particles:
        # Calculate time for collision with vertical (left/right) walls
        t_x, wall_x = calculate_wall_collision(p, L_x, radius, p.xvel, False)
        # Calculate time for collision with horizontal (lower/upper) walls
        t_y, wall_y = calculate_wall_collision(p, L_y, radius, p.yvel, True)
        
        # Add the earlier wall collision to the array
        if t_x < t_y:
            collision_array.append(Collision(t_x, wall_x, p.id, "wall"))
        else:
            collision_array.append(Collision(t_y, wall_y, p.id, "wall"))

    # Calculate potential particle-particle collisions and set them to MAX initially
    for p1, p2 in combinations(range(n_particles), 2):
        collision_array.append(Collision(MAX_TIME, p1, p2, "particle"))

    return collision_array

# Initialize the collision array
collision_array = initialize_collision_array(particles, L_x, L_y, radius)

# Display the initial collision array for verification
for collision in collision_array:
    time = f"{collision.time:.6f}" if collision.time != MAX_TIME else "MAX"
    print(f"{time}: {collision.part1}, {collision.part2}, ({collision.ctype})")



# Function to sort the collision array by time
def sort_collision_array(collision_array):
    # Insertion sort is used here for simplicity
    for i in range(1, len(collision_array)):
        key = collision_array[i]
        j = i - 1
        # Move elements of collision_array[0..i-1], that are greater than key,
        # to one position ahead of their current position
        while j >= 0 and key.time < collision_array[j].time:
            collision_array[j + 1] = collision_array[j]
            j -= 1
        collision_array[j + 1] = key

# Sort the initial collision array
sort_collision_array(collision_array)

# Display the sorted initial collision array for verification
for collision in collision_array:
    time = f"{collision.time:.6f}" if collision.time != MAX_TIME else "MAX"
    print(f"{time}: {collision.part1}, {collision.part2}, ({collision.ctype})")


def process_wall_collision(particle, collision_time, wall_id, L_x, L_y):
    # Update the particle's position to the time of collision
    particle.xpos += particle.xvel * (collision_time - particle.last_update)
    particle.ypos += particle.yvel * (collision_time - particle.last_update)
    particle.last_update = collision_time  # Update the last update time to the collision time

    # If the collision is with a vertical wall, reverse the x-velocity
    if wall_id == LEFT_WALL or wall_id == RIGHT_WALL:
        particle.xpos = L_x - radius if particle.xvel > 0 else radius  # Correct position to wall position
        particle.xvel = -particle.xvel
        particle.cwall += 1  # Increment the wall collision count

    # If the collision is with a horizontal wall, reverse the y-velocity
    if wall_id == LOWER_WALL or wall_id == UPPER_WALL:
        particle.ypos = L_y - radius if particle.yvel > 0 else radius  # Correct position to wall position
        particle.yvel = -particle.yvel
        particle.cwall += 1  # Increment the wall collision count


# Function to calculate time to collision with a wall
def time_to_wall(particle, wall_id, L_x, L_y):
    if wall_id == LEFT_WALL:  # Collision with left wall
        if particle.xvel < 0:
            return (radius - particle.xpos) / particle.xvel
    elif wall_id == RIGHT_WALL:  # Collision with right wall
        if particle.xvel > 0:
            return (L_x - radius - particle.xpos) / particle.xvel
    elif wall_id == LOWER_WALL:  # Collision with lower wall
        if particle.yvel < 0:
            return (radius - particle.ypos) / particle.yvel
    elif wall_id == UPPER_WALL:  # Collision with upper wall
        if particle.yvel > 0:
            return (L_y - radius - particle.ypos) / particle.yvel
    return MAX_TIME

# Function to calculate time to collision between two particles
def time_to_particle(p1, p2):
    # Calculate relative position and velocity
    dx = p2.xpos - p1.xpos
    dy = p2.ypos - p1.ypos
    dvx = p2.xvel - p1.xvel
    dvy = p2.yvel - p1.yvel

    # Calculate discriminant
    a = dvx**2 + dvy**2
    b = 2 * (dx * dvx + dy * dvy)
    c = dx**2 + dy**2 - (2*radius)**2
    discriminant = b**2 - 4 * a * c

    # If discriminant is negative, no collision is possible
    if discriminant < 0:
        return MAX_TIME
    
    # Calculate time to collision
    t_collision = (-b - np.sqrt(discriminant)) / (2 * a)
    if t_collision < 0:
        return MAX_TIME  # Collision is in the past, not possible
    return t_collision




It seems there is a discrepancy in the calculation. The expected next collision should be between particles 1 and 2 at approximately time 3.374136 seconds. However, the collision array provided shows a different set of nearest collision times.

Let's verify the calculations once more, particularly focusing on the collision between particles 1 and 2, to ensure we have the correct nearest collision time. I'll go ahead and recompute the times for potential collisions between particles 1 and 2 based on their current state and velocities.

Finished working
Hide work

python
Copy code
# Function to compute the time of the next collision between two particles
def compute_particle_collision_time(p1, p2, radius):
    # Relative position
    dx = p1['xpos'] - p2['xpos']
    dy = p1['ypos'] - p2['ypos']
    
    # Relative velocity
    dvx = p1['xvel'] - p2['xvel']
    dvy = p1['yvel'] - p2['yvel']
    
    # Coefficients of the quadratic equation (a*t^2 + b*t + c = 0)
    a = dvx**2 + dvy**2
    b = 2 * (dx*dvx + dy*dvy)
    c = dx**2 + dy**2 - (2*radius)**2
    
    # Discriminant
    discriminant = b**2 - 4*a*c
    
    if discriminant < 0:
        # Particles are not going to collide
        return None
    
    # Time to collision
    discriminant_sqrt = np.sqrt(discriminant)
    t1 = (-b - discriminant_sqrt) / (2*a)
    t2 = (-b + discriminant_sqrt) / (2*a)
    
    # We take the smallest positive time
    if t1 > 0 and t2 > 0:
        return min(t1, t2)
    elif t1 > 0:
        return t1
    elif t2 > 0:
        return t2
    else:
        # Collision is in the past
        return None