// Find diameter using All-pairs shortest paths
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "graph.h"

// #define DEBUG

//===================== Function Prototypes ===========================
void PrepareGraph(int* mygraph);
int calcIterations();
void CalcShortestPaths(int* mygraph);
int CalcDiameter(int* mygraph);
void PrintGraph(int* mygraph);

//========================= Main Function ===========================
int main() {
// initialize the graph here; note it is a 1-dimensional array	
	int* mygraph = (int *) malloc(sizeof(int) * NVERT * NVERT);
	getGraph(mygraph);

	int i, j, k;

#ifdef DEBUG
	printf("\nThe adjacency matrix is\n");
	PrintGraph(mygraph);
#endif

// prepare the graph/array for the all-pairs algorithm
	PrepareGraph(mygraph);

#ifdef DEBUG
	printf("\nThe matrix for calculating all-pairs shortest paths is\n");
	PrintGraph(mygraph);
#endif
	
// set up outer iteration
	int niterations = calcIterations();
#ifdef DEBUG
	printf("number of iterations = %d\n",niterations);
#endif

// calculate shortest paths	
	for(int iiteration = 0; iiteration < niterations; iiteration++) {
#ifdef DEBUG
		printf("\nouter iteration %d:\n",iiteration);
#endif
		CalcShortestPaths(mygraph);
	} 
	
#ifdef DEBUG
	printf("\nThe graph diameter is ");
#endif
	printf("%d\n",CalcDiameter(mygraph));

	free(mygraph);
	return 0;
}

//===================== Function Definitions ===========================

// store "Infinity" values (INT_MAX) where there are no edges, other than from each edge to itself
void PrepareGraph(int* mygraph) {
}

int calcIterations() {
	return 0;
}

// calculate the lengths of all shortest paths using the "matrix multiplication" method
void CalcShortestPaths(int* mygraph) {
}

// calculate the diameter of the graph as the maximum of the minimum distances
int CalcDiameter(int* mygraph) {
	return 0;
}

// print all array entries as an adjacency matrix
void PrintGraph(int* mygraph) {
	for (int i=0; i<NVERT; i++) {
		for (int j=0; j<NVERT; j++) {
			printf("%d ", mygraph[i*NVERT+j]);
		}
		printf("\n");
	}
}