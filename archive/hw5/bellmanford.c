#include "bellmanford.h"

#include <math.h>
#include <stdbool.h>

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
        // Output shortest path for a specified destination
        printf("%d: %.5f; ", destination, distances[destination]);
        // Backtrack from the destination to source to print the path
        int current = destination;
        while (current != -1) {
            printf("%d ", current);
            current = predecessors[current];
        }
        printf("\n");
    } else {
        // Output shortest paths for all vertices
        for (int i = 0; i < numvertices; i++) {
            printf("%d: %.5f; ", i, distances[i]);
            int current = i;
            while (current != -1) {
                printf("%d ", current);
                current = predecessors[current];
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
