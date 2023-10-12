#include <stdio.h>
#include <stdlib.h>
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

// Function declarations
Graph* createGraph(int vertices);
Node* createNode(int vertex, int label);
void addEdge(Graph* graph, int src, int dest, int label);
void DFS(Graph* graph, int vertex, int maxLabel, int* visited);
void freeGraph(Graph* graph);
Graph* read_graph_from_file(const char* filename);
int count_components(Graph* graph);

int main(int argc, char* argv[]) {
  if (argc != 3) {
    printf("Usage: ./program_name <filename> <max_label>\n");
    return 1;
  }

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

  Graph* graph = read_graph_from_file(filename);
  if (!graph) { return; }

  int components = count_components(graph);
  printf("%d\n", components);
  free_graph(graph)
  return 0;
}

Graph* read_graph_from_file(const char* filename) {
  FILE* file = fopen(filename, "r");
  if (!file) {
    printf("Error");
    return NULL;
  }
  int numVertices;
  fscanf(file, "%d", &numVertices);
  Graph* graph = malloc(sizeof(Graph));
  if (!graph) {
    fclose(file);
    return NULL;
  }
  graph->numVertices = numVertices;
  graph->adjLists = malloc(sizeof(numVertices * Node*));
  for (int i = 0; i < numVertices; i++) {
    graph->adjLists[i] = NULL;
  }
  int src, dest, label;
  while (fscanf(file, "%d %d %d", &src, &dest, &label) == 3) {
    if (label <= maxLabel) {
      add_edge(graph, src, dest, label);
    }
  }
  fclose(file);
  return graph;
}

void add_edge(Graph* graph, int src, int dest, int label) {
  Node* newNode = malloc(sizeof(Node));
  newNode->vertex = dest;
  newNode->label = label;
  newNode->next = graph->adjList[src];
  graph->adjList[src] = newNode;

  newNode = malloc(sizeof(Node));
  newNode->vertex = src;
  newNode->label = label;
  newNode->next = graph->adjList[dest];
  graph->adjList[dest] = newNode;
}  

void free_graph(Graph* graph) {
  if (!graph) return;
  for (int i = 0; i < numVertices; i++) {
    Node* temp = graph->adjLists[i];
    Node* nextNode = graph->adjLists[i]->next;
    while (temp->next) {
      temp->next = nextNode->next;
      free(nextNode);
      nextNode = temp->next;
    }
    free(temp);
    
//     while (temp) { 
//       Node* nextNode = temp->next;
//       free(temp);
//       temp = nextNode;
//     }
//   }
  free(graph->adjList);
  free(graph);
}

void dfs(Graph* graph, int vertex, bool visited[], int maxLabel) {
  visited[vertex] = true;
  Node* temp = graph->adjLists[vertex];
  while (temp) {
    if (temp->label <= maxLabel && !is_visited(visited, temp->vertex)) {
      dfs(graph, temp->vertex, visited, maxLabel);
    }
    temp = temp->next;
  }
}
bool is_visited(bool visited[], int vertex) {
  return visited[vertex];
}

int count_components(Graph* graph) {
  int count = 0;
  int numVertices = graph->numVertices;
  bool* visited = malloc(numVertices * sizeof(bool));
  for (int i = 0; i < numVertices; i++) { visited[i] = false; }
  for (int v = 0; v < numVertices; v++) {
    if (!visited[v]) {
      dfs(graph, v, visited);
      count++;
    }
  }
  free(visited);
  return count;
}

// Graph* createGraph(int numVertices) {
//   Graph* graph = malloc(sizeof(Graph));
//   graph->numVertices = numVertices;
//   graph->adjLists = malloc(numVertices * sizeof(Node*));
//   for (int i = 0; i < numVertices; i++) { graph->adjLists[i] = NULL; }
//   return graph;
// }

// Node* createNode(int vertex, int label) {
//   Node* newNode = malloc(sizeof(Node));
//   newNode->vertex = vertex;
//   newNode->label = label;
//   newNode->next = NULL;
//   return newNode;
// }

// void addEdge(Graph* graph, int src, int dest, int label) {
//   Node* newNode = createNode(dest, label);
//   newNode->next = graph->adjLists[src];
//   graph->adjLists[src] = newNode;
//   newNode = createNode(src, label);
//   newNode->next = graph->adjLists[dest];
//   graph->adjLists[dest] = newNode;
// }

// void DFS(Graph* graph, int vertex, int* visited, int maxLabel) {
//   Node* temp = graph->adjLists[vertex];
//   visited[vertex] = 1;
//   while (temp) {
//     if (!visited[temp->vertex] && temp->label <= maxLabel) {
//       DFS(graph, temp->vertex, visited, maxLabel);
//     }
//     temp = temp->next;
//   }
// }

// int countConnectedComponents(Graph* graph, int maxLabel) {
//   int* visited = calloc(graph->numVertices, sizeof(int));
//   int components = 0;
//   for (int v = 0; v < graph->numVertices; v++) {
//     if (!visited[v]) {
//       DFS(graph, v, visited, maxLabel);
//       components++;
//     }
//   }
//   free(visited);
//   return components;
// }

// void freeGraph(Graph* graph) {
//   for (int i = 0; i < graph->numVertices; i++) {
//     Node* temp = graph->adjLists[i];
//     while (temp) {
//       Node* next = temp->next;
//       free(temp);
//       temp = next;
//     }
//   }
//   free(graph->adjLists);
//   free(graph);
// }
