#include "bellmanford.h"

void readGraph(const char* filename, int* numvertices, int* numedges, edge** edges) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    fscanf(file, "%d %d", numvertices, numedges);
    *edges = malloc(*numedges * sizeof(edge));

    for (int i = 0; i < *numedges; i++) {
        fscanf(file, "%d %d %lf", &(*edges)[i].from, &(*edges)[i].to, &(*edges)[i].weight);
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
    bool hasChanged;

    for (int i = 0; i < numvertices - 1; i++) {
        hasChanged = false;

        #pragma omp parallel for
        for (int j = 0; j < numedges; j++) {
            int u = edges[j].from;
            int v = edges[j].to;
            double weight = edges[j].weight;

            if (distances[u] + weight < distances[v]) {
                distances[v] = distances[u] + weight;
                predecessors[v] = u;
                hasChanged = true;
            }
        }

        if (!hasChanged) {
            break;
        }
    }
}

void outputResult(int numvertices, int destination, double* distances, int* predecessors) {
    if (destination >= 0 && destination < numvertices) {
        // Temporary array to store the path
        int *path = (int *)malloc(numvertices * sizeof(int));
        int pathLength = 0;

        // Output shortest path for a specified destination
        printf("%d: %.5f; ", destination, distances[destination]);

        // Backtrack from the destination to source and store the path
        int current = destination;
        while (current != -1) {
            path[pathLength++] = current;
            current = predecessors[current];
        }

        // Reverse the order and print the path
        for (int i = pathLength - 1; i >= 0; i--) {
            printf("%d ", path[i]);
        }
        printf("\n");

        // Free the allocated memory
        free(path);
    } else {
        // Output shortest paths for all vertices
        for (int i = 0; i < numvertices; i++) {
            // Temporary array to store the path
            int *path = (int *)malloc(numvertices * sizeof(int));
            int pathLength = 0;

            printf("%d: %.5f; ", i, distances[i]);

            // Backtrack from the destination to source and store the path
            int current = i;
            while (current != -1) {
                path[pathLength++] = current;
                current = predecessors[current];
            }

            // Reverse the order and print the path
            for (int j = pathLength - 1; j >= 0; j--) {
                printf("%d ", path[j]);
            }
            printf("\n");

            // Free the allocated memory
            free(path);
        }
    }
}

void cleanup(edge* edges, double* distances, int* predecessors) {
    free(edges);
    free(distances);
    free(predecessors);
}
