import sys
import time

def read_graph(filename):
    with open(filename, 'r') as file:
        num_vertices, num_edges = map(int, file.readline().split())
        edges = [list(map(float, line.split())) for line in file]
    return num_vertices, edges

def bellman_ford(num_vertices, edges, source=0):
    distances = [float('inf')] * num_vertices
    predecessors = [-1] * num_vertices
    distances[source] = 0

    for _ in range(num_vertices - 1):
        for u, v, weight in edges:
            if distances[int(u)] + weight < distances[int(v)]:
                distances[int(v)] = distances[int(u)] + weight
                predecessors[int(v)] = int(u)

    for u, v, weight in edges:
        if distances[int(u)] + weight < distances[int(v)]:
            raise ValueError("Graph contains a negative weight cycle")

    return distances, predecessors

def format_path(predecessors, destination):
    path = []
    while destination != -1:
        path.append(destination)
        destination = predecessors[destination]
    return path[::-1]

def output_results(num_vertices, distances, predecessors, destination=None):
    if destination is not None:
        destinations = [destination]
    else:
        destinations = range(num_vertices)

    for vertex in destinations:
        path = format_path(predecessors, vertex)
        path_str = ' '.join(map(str, path))
        print(f"{vertex}: {distances[vertex]:.5f}; {path_str}")

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 bellmanford.py <filename> <num_threads> [destination_vertex]")
        sys.exit(1)

    filename = sys.argv[1]
    # num_threads = int(sys.argv[2])  # Not used in Python, but kept for consistency with command line arguments
    destination = int(sys.argv[3]) if len(sys.argv) > 3 else None

    num_vertices, edges = read_graph(filename)
    distances, predecessors = bellman_ford(num_vertices, edges)
    output_results(num_vertices, distances, predecessors, destination)

if __name__ == "__main__":
    main()
