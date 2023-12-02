#ifndef CONNECTION_H
#define CONNECTION_H

#include <stdbool.h>

// Node struct for representing individual nodes in the adjacency list.
typedef struct Node {
    int vertex;         // The vertex number/index
    int label;          // Edge label associated with the vertex
    struct Node* next;  // Pointer to the next node in the adjacency list
} Node;

// Graph struct for representing the entire graph.
typedef struct Graph {
    int numVertices;    // Total number of vertices in the graph
    Node** adjLists;    // Array of pointers, each pointing to an adjacency list for a vertex
} Graph;

// Function to create a new graph with a given number of vertices.
Graph* createGraph(int numVertices);

// Function to create a new node with the specified vertex number and label.
Node* createNode(int vertex, int label);

// Function to add an edge to the graph between the given source and destination vertices.
void addEdge(Graph* graph, int src, int dest, int label);

// Function to deallocate memory used by the graph.
void freeGraph(Graph* graph);

// Function for Depth-First Search (DFS) starting from a given vertex.
void dfs(Graph* graph, int vertex, bool visited[], int maxLabel);

// Function to check if a given vertex has been visited during DFS.
bool isVisited(bool visited[], int vertex);

// Function to count the number of connected components in the graph.
int countComponents(Graph* graph, int maxLabel);

// Utility function to read a graph from a file.
Graph* readGraphFromFile(const char* filename, int maxLabel);

#endif // CONNECTION_H
