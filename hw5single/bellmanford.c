#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>

// Structure definition
typedef struct {
    int from, to;
    double weight;
} edge;

void readGraph(const char* filename, int* numvertices, int* numedges, edge** edges) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    if (fscanf(file, "%d %d", numvertices, numedges) != 2) {
        perror("Error reading graph size");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    *edges = malloc(*numedges * sizeof(edge));
    if (*edges == NULL) {
        perror("Failed to allocate memory for edges");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < *numedges; i++) {
        if (fscanf(file, "%d %d %lf", &(*edges)[i].from, &(*edges)[i].to, &(*edges)[i].weight) != 3) {
            perror("Error reading edge data");
            fclose(file);
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);
}

void initialize(int numvertices, double* distances, int* predecessors) {
    for (int i = 0; i < numvertices; i++) {
        distances[i] = (i == 0) ? 0 : INFINITY;
        predecessors[i] = -1;
    }
}

void bellmanFord(int numvertices, int numedges, edge* edges, double* distances, int* predecessors) {
    bool globalChange;

    for (int i = 0; i < numvertices - 1; i++) {
        bool localChange = false;

        #pragma omp parallel for reduction(|:localChange)
        for (int j = 0; j < numedges; j++) {
            int u = edges[j].from;
            int v = edges[j].to;
            double weight = edges[j].weight;

            if (distances[u] + weight < distances[v]) {
                distances[v] = distances[u] + weight;
                predecessors[v] = u;
                localChange = true;
            }
        }

        globalChange |= localChange;
        if (!globalChange) {
            break;
        }
    }
}

void outputResult(int numvertices, int destination, double* distances, int* predecessors) {
    int path[numvertices];
    int pathLength = 0;

    if (destination >= 0 && destination < numvertices) {
        // Output shortest path for a specified destination
        printf("%d: %.5f; ", destination, distances[destination]);
        
        int current = destination;
        while (current != -1) {
            path[pathLength++] = current;
            current = predecessors[current];
        }

        for (int i = 0; i < pathLength; i++) {
            printf("%d ", path[i]);
        }
        printf("\n");

    } else {
        // Output shortest paths for all vertices
        for (int i = 0; i < numvertices; i++) {
            printf("%d: %.5f; ", i, distances[i]);
            
            int current = i;
            pathLength = 0;
            while (current != -1) {
                path[pathLength++] = current;
                current = predecessors[current];
            }

            for (int j = 0; j < pathLength; j++) {
                printf("%d ", path[j]);
            }
            printf("\n");
        }
    }
}

void cleanup(edge* edges, double* distances, int* predecessors) {
    free(edges);
    free(distances);
    free(predecessors);
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <filename> <num_threads> [destination_vertex]\n", argv[0]);
        return 1;
    }

    const char* filename = argv[1];
    int numThreads = atoi(argv[2]);
    if (numThreads <= 0) {
        fprintf(stderr, "Number of threads must be positive\n");
        return 1;
    }

    int destination = (argc > 3) ? atoi(argv[3]) : -1;

    int numvertices, numedges;
    edge* edges = NULL;

    readGraph(filename, &numvertices, &numedges, &edges);
    double* distances = malloc(numvertices * sizeof(double));
    if (distances == NULL) {
        perror("Failed to allocate memory for distances");
        exit(EXIT_FAILURE);
    }
    int* predecessors = malloc(numvertices * sizeof(int));
    if (predecessors == NULL) {
        perror("Failed to allocate memory for predecessors");
        free(distances);
        exit(EXIT_FAILURE);
    }

    initialize(numvertices, distances, predecessors);
    omp_set_num_threads(numThreads);

    double startTime = omp_get_wtime();
    bellmanFord(numvertices, numedges, edges, distances, predecessors);
    double endTime = omp_get_wtime();

    outputResult(numvertices, destination, distances, predecessors);

    printf("Elapsed time: %.5f seconds\n", endTime - startTime);

    cleanup(edges, distances, predecessors);

    return 0;
}
 
