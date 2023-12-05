#ifndef BELLMANFORD_H
#define BELLMANFORD_H

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

// Function declarations
void readGraph(const char* filename, int* numvertices, int* numedges, edge** edges);
void initialize(int numvertices, double* distances, int* predecessors);
void bellmanFord(int numvertices, int numedges, edge* edges, double* distances, int* predecessors);
void outputResult(int numvertices, int destination, double* distances, int* predecessors);
void cleanup(edge* edges, double* distances, int* predecessors);

#endif // BELLMANFORD_H
