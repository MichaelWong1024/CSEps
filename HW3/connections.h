#ifndef CONNECTION_H
#define CONNECTION_H

#include <stdbool.h>

typedef struct Node {
    int vertex;
    int label;
    struct Node* next;
} Node;

typedef struct Graph {
    int numVertices;
    Node** adjLists;
} Graph;

// Function declarations for graph operations
Graph* createGraph(int numVertices);
Node* createNode(int vertex, int label);
void addEdge(Graph* graph, int src, int dest, int label);
void freeGraph(Graph* graph);

// DFS related functions
void dfs(Graph* graph, int vertex, bool visited[], int maxLabel);
bool isVisited(bool visited[], int vertex);
int countComponents(Graph* graph, int maxLabel);

// Utility functions
Graph* readGraphFromFile(const char* filename, int maxLabel);

#endif // CONNECTION_H
