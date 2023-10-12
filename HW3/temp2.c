#include <stdio.h>
#include <stdlib.h>

typedef struct Node {
    int vertex;
    int label;
    struct Node* next;
} Node;

typedef struct Graph {
    int numVertices;
    Node** adjLists;
} Graph;

// Function declarations
Graph* createGraph(int vertices);
Node* createNode(int vertex, int label);
void addEdge(Graph* graph, int src, int dest, int label);
void DFS(Graph* graph, int vertex, int maxLabel, int* visited);
void freeGraph(Graph* graph);

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("Usage: ./program_name <filename> <max_label>\n");
        return 1;
    }

    // Read input file and max label
    char* filename = argv[1];
    int maxLabel = atoi(argv[2]);
    if (maxLabel != 1 && maxLabel != 2) {
        printf("Invalid max label. Defaulting to 2.\n");
        maxLabel = 2;
    }

    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error in opening the file.\n");
        return 1;
    }

    int n, u, v, label;
    fscanf(file, "%d", &n);

    Graph* graph = createGraph(n);
    while (fscanf(file, "%d %d %d", &u, &v, &label) != EOF) {
        addEdge(graph, u, v, label);
        addEdge(graph, v, u, label);
    }

    fclose(file);

    int count = 0;
    int* visited = (int*)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        visited[i] = 0;
    }

    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            DFS(graph, i, maxLabel, visited);
            count++;
        }
    }

    printf("%d\n", count);

    free(visited);
    freeGraph(graph);
    return 0;
}

Graph* createGraph(int vertices) {
    Graph* graph = malloc(sizeof(Graph));
    graph->numVertices = vertices;
    graph->adjLists = malloc(vertices * sizeof(Node*));
    for (int i = 0; i < vertices; i++) {
        graph->adjLists[i] = NULL;
    }
    return graph;
}

Node* createNode(int vertex, int label) {
    Node* newNode = malloc(sizeof(Node));
    newNode->vertex = vertex;
    newNode->label = label;
    newNode->next = NULL;
    return newNode;
}

void addEdge(Graph* graph, int src, int dest, int label) {
    Node* newNode = createNode(dest, label);
    newNode->next = graph->adjLists[src];
    graph->adjLists[src] = newNode;
}

void DFS(Graph* graph, int vertex, int maxLabel, int* visited) {
    Node* temp = graph->adjLists[vertex];
    visited[vertex] = 1;
    while (temp) {
        int connectedVertex = temp->vertex;
        if (!visited[connectedVertex] && temp->label <= maxLabel) {
            DFS(graph, connectedVertex, maxLabel, visited);
        }
        temp = temp->next;
    }
}

void freeGraph(Graph* graph) {
    for (int i = 0; i < graph->numVertices; i++) {
        Node* temp = graph->adjLists[i];
        while (temp) {
            Node* next = temp->next;
            free(temp);
            temp = next;
        }
    }
    free(graph->adjLists);
    free(graph);
}
