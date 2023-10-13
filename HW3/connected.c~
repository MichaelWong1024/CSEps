#include <stdio.h>
#include <stdlib.h>
#include "connections.h"

Graph* readGraphFromFile(const char* filename, int maxLabel) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening the file: %s\n", filename);
        return NULL;
    }

    int numVertices;
    fscanf(file, "%d", &numVertices);
    
    Graph* graph = createGraph(numVertices);
    if (!graph) {
        printf("Error allocating memory for the graph.\n");
        fclose(file);
        return NULL;
    }

    int src, dest, label;
    while (fscanf(file, "%d %d %d", &src, &dest, &label) == 3) {
        if (label <= maxLabel) {
            addEdge(graph, src, dest, label);
        }
    }
    fclose(file);
    return graph;
}

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

Node* createNode(int vertex, int label) {
    Node* newNode = malloc(sizeof(Node));
    if (!newNode) return NULL;

    newNode->vertex = vertex;
    newNode->label = label;
    newNode->next = NULL;
    return newNode;
}

void addEdge(Graph* graph, int src, int dest, int label) {
    Node* newNode = createNode(dest, label);
    newNode->next = graph->adjLists[src];
    graph->adjLists[src] = newNode;

    newNode = createNode(src, label);
    newNode->next = graph->adjLists[dest];
    graph->adjLists[dest] = newNode;
}  

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

bool isVisited(bool visited[], int vertex) {
    return visited[vertex];
}

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
