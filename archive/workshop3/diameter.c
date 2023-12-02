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
  for (int i = 0; i < NVERT; i++) {
    for (int j = 0; j < NVERT; j++) {
      if (i != j && mygraph[i * NVERT + j] == 0) {
	mygraph[i * NVERT + j] = INT_MAX;
      }
    }
  }
}

int calcIterations() {
  double iterations = ceil(log(NVERT - 1) / log(2));
  return iterations;
}

// calculate the lengths of all shortest paths using the "matrix multiplication" method
void CalcShortestPaths(int* mygraph) {
  int* temp = (int *) malloc(sizeof(int) * NVERT * NVERT); // Create a temporary matrix
  int index_ij, index_ik, index_kj; // Indices for the elements of the 1D representation of the matrix

  for (int i = 0; i < NVERT; i++) {
    for (int j = 0; j < NVERT; j++) {
      index_ij = i * NVERT + j;
      temp[index_ij] = mygraph[index_ij];

      for (int k = 0; k < NVERT; k++) {
	index_ik = i * NVERT + k;
	index_kj = k * NVERT + j;

	if (mygraph[index_ik] != INT_MAX && mygraph[index_kj] != INT_MAX) {
	  if (temp[index_ij] > (long)mygraph[index_ik] + mygraph[index_kj]) {
	    temp[index_ij] = mygraph[index_ik] + mygraph[index_kj];
	  }
	}
      }
    }
  }
  for (int i = 0; i < NVERT * NVERT; i++) {
    mygraph[i] = temp[i];
  }
  free(temp);
}

// calculate the diameter of the graph as the maximum of the minimum distances
int CalcDiameter(int* mygraph) {
  int diameter = 0;
  for (int i = 0; i < NVERT; i++) {
    for (int j = 0; j < NVERT; j++) {
      if (mygraph[i * NVERT + j] != INT_MAX && mygraph[i * NVERT +j] > diameter) {
	diameter = mygraph[i * NVERT + j];
      }
    }
  }
  return diameter;
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

      
