#include <stdio.h>
#include <stdlib.h>
#include "connected.h"

// Read the graph from the provided file and return a graph representation.
// Only edges with labels less than or equal to maxLabel are considered.
Graph* readGraphFromFile(const char* filename, int maxLabel) {
    // Open the provided file for reading.
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening the file: %s\n", filename);
        return NULL;
    }

    int numVertices;
    // Read the number of vertices from the file.
    fscanf(file, "%d", &numVertices);

    // Create a new graph structure.
    Graph* graph = createGraph(numVertices);
    if (!graph) {
        printf("Error allocating memory for the graph.\n");
        fclose(file);
        return NULL;
    }

    int src, dest, label;
    // Read each edge from the file and add it to the graph.
    while (fscanf(file, "%d %d %d", &src, &dest, &label) == 3) {
        if (label <= maxLabel) {
            addEdge(graph, src, dest, label);
        }
    }

    fclose(file);
    return graph;
}

// Create a new graph structure with the specified number of vertices.
Graph* createGraph(int numVertices) {
    Graph* graph = malloc(sizeof(Graph));
    if (!graph) return NULL;

    graph->numVertices = numVertices;
    graph->adjLists = malloc(numVertices * sizeof(Node*));

    if (!graph->adjLists) {
        free(graph);
        return NULL;
    }

    for (int i = 0; i < numVertices; i++) {
        graph->adjLists[i] = NULL;
    }
    return graph;
}

// Create a new node for a vertex with the given label.
Node* createNode(int vertex, int label) {
    Node* newNode = malloc(sizeof(Node));
    if (!newNode) return NULL;

    newNode->vertex = vertex;
    newNode->label = label;
    newNode->next = NULL;
    return newNode;
}

// Add an edge between the given source and destination vertices with the specified label.
void addEdge(Graph* graph, int src, int dest, int label) {
    Node* newNode = createNode(dest, label);
    newNode->next = graph->adjLists[src];
    graph->adjLists[src] = newNode;

    // Since the graph is undirected, we also add a reverse edge.
    newNode = createNode(src, label);
    newNode->next = graph->adjLists[dest];
    graph->adjLists[dest] = newNode;
}

// Deallocate memory used by the graph.
void freeGraph(Graph* graph) {
    for (int i = 0; i < graph->numVertices; i++) {
        Node* current = graph->adjLists[i];
        while (current) {
            Node* temp = current;
            current = current->next;
            free(temp);
        }
    }
    free(graph->adjLists);
    free(graph);
}

// Depth-first search, starting from the given vertex.
void dfs(Graph* graph, int vertex, bool visited[], int maxLabel) {
    visited[vertex] = true;
    Node* temp = graph->adjLists[vertex];
    while (temp) {
        if (temp->label <= maxLabel && !isVisited(visited, temp->vertex)) {
            dfs(graph, temp->vertex, visited, maxLabel);
        }
        temp = temp->next;
    }
}

// Check if a vertex has been visited during DFS.
bool isVisited(bool visited[], int vertex) {
    return visited[vertex];
}

// Count the number of connected components in the graph.
int countComponents(Graph* graph, int maxLabel) {
    int count = 0;
    bool* visited = calloc(graph->numVertices, sizeof(bool));
    if (!visited) return -1;

    for (int v = 0; v < graph->numVertices; v++) {
        if (!visited[v]) {
            dfs(graph, v, visited, maxLabel);
            count++;
        }
    }
    free(visited);
    return count;
}
