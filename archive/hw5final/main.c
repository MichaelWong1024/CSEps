#include "bellmanford.h"

int main(int argc, char* argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <filename> <num_threads> [destination_vertex]\n", argv[0]);
        return 1;
    }

    const char* filename = argv[1];
    int numThreads = atoi(argv[2]);
    int destination = (argc > 3) ? atoi(argv[3]) : -1;

    int numvertices, numedges;
    edge* edges = NULL;

    readGraph(filename, &numvertices, &numedges, &edges);
    double* distances = malloc(numvertices * sizeof(double));
    int* predecessors = malloc(numvertices * sizeof(int));
    initialize(numvertices, distances, predecessors);
    omp_set_num_threads(numThreads);

    double start_time = omp_get_wtime();  // Start time measurement
    bellmanFord(numvertices, numedges, edges, distances, predecessors);
    double end_time = omp_get_wtime();  // End time measurement

    outputResult(numvertices, destination, distances, predecessors);
    printf("%.5f\n", end_time - start_time);  // Print elapsed time

    cleanup(edges, distances, predecessors);
}