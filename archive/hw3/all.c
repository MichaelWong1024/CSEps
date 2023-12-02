//
//  all.c
//  6010hw3
//
//  Created by Zhikai Xu on 2023/10/4.
//

#include "all.h"
#include <stdio.h>
#include <stdlib.h>

// initialize a graph using the given number of verices
Graph* createGraph(int numVertices) {
    Graph* newGraph = malloc(sizeof(Graph));
    newGraph->numVertices = numVertices;
    // initialize the adjacency list
    newGraph->adjlist = (Node**) malloc(numVertices * sizeof(Node*));
    for (int i = 0; i < numVertices; i++) {
        newGraph->adjlist[i] = NULL;
    }
    // initialize the the from-to 1D array
    // here a adjacency matrix (squeezed to be 1D array) is used to represent the from-to relationship
    newGraph->fromTo = (int*) malloc(numVertices * numVertices * sizeof(int));
    for (int i = 0; i < numVertices * numVertices; i++) {
        newGraph->fromTo[i] = 0;
    }
    // the max edge label given by the user
    newGraph->maxLabel = -1;
    return newGraph;
}

// initialize the node of the adjacency list
Node* createNode(int vertex) {
    Node* newNode = malloc(sizeof(Node));
    newNode->vertex = vertex;
    newNode->next = NULL;
    return newNode;
}

// add a edge into the adjacency list and the from-to array
void addEdge(Graph* graph, int from, int vertex, int edgeLabel) {
    Node* newNode = createNode(vertex);
    Node* temp = (graph->adjlist)[from];
    newNode->next = temp;
    // add the node into the adjacency list
    graph->adjlist[from] = newNode;
    // add the from-to relationship into from-to array
    graph->fromTo[from * graph->numVertices + vertex] = edgeLabel;
}

// read the graph file and store the graph in a 1D-array
Graph* readGraphFromFile(char* fileName, int edgeLabelGiven) {
    if (edgeLabelGiven != 1 && edgeLabelGiven != 2) {
        edgeLabelGiven = 2;
    }
    FILE* fp;
    int numV;
    fp = fopen(fileName, "r");
    
    // read the number of vertices
    fscanf(fp, "%d", &numV);
    
    // allocate memory for the graph to create the graph
    Graph* newGraph = createGraph(numV);
    int data;
    int count = 0; // use to tell the meaning of the data
    int from = -1;
    int to = -1;
    int edgeLabel = -1;
    while (fscanf(fp, "%d", &data) != EOF) {
        if (count % 3 == 0) {
            from = data; // read the from vertex
        } else if (count % 3 == 1) {
            to = data; // read the to vertex
        } else if (count % 3 == 2){
            edgeLabel = data; // read the edge label
        } else {
            printf("Error in readGraphFromFile()!");
            break;
        }
        // after read a line, add the edge into adjacency list and from-to array
        if (count % 3 == 2 && edgeLabel <= edgeLabelGiven) {
            addEdge(newGraph, from, to, edgeLabel);
        }
        count++;
    }
    fclose(fp);
    return newGraph;
}

// initialize a stack for the depth-first search
Stack* createStack(int capacity) {
    Stack* stack = malloc(sizeof(Stack));
    stack->capacity = capacity;
    stack->top = -1;
    stack->array = malloc(capacity * sizeof(int));
    return stack;
}

// return 1 if stack is full, otherwise return 0
int isFull(Stack* stack) {
    return stack->top == stack->capacity - 1;
}

// return 1 if stack is empty, otherwise return 0
int isEmpty(Stack* stack) {
    return stack->top == -1;
}

// push a item onto the top of the stack
void push(Stack* stack, int item) {
    if (isFull(stack)) {
        return;
    }
    stack->top++;
    (stack->array)[stack->top] = item;
}

// remove a item from the top of the stack, and return the removed item value
int pop(Stack* stack) {
    if (isEmpty(stack)) {
        return -1;
    }
    return stack->array[stack->top--];
}

// using dfs to count the connected components
int dfs(Graph* graph) {
    int count = 0;
    int* visited = malloc(graph->numVertices * sizeof(int));
    for (int i = 0; i < graph->numVertices; i++) {
        visited[i] = 0;
    }
    Stack* stack = createStack(graph->numVertices);
    for (int i = 0; i < graph->numVertices; i++) {
        // traverse it using dfs
        if (visited[i] == 0) {
            push(stack, i);
            visited[i] = 1;
            while (!isEmpty(stack)) {
                int from = pop(stack);
                Node* to = graph->adjlist[from];
                // iterate through each linked list of each vertex
                while (to != NULL) {
                    if (visited[to->vertex] == 0) {
                        push(stack, to->vertex);
                        visited[to->vertex] = 1;
                    }
                    to = to->next;
                }
            }
            // every time dfs finishs, we find a new connected component
            count++;
        }
    }
    
    // free all the allocated memory
    free(visited);
    free(stack->array);
    free(stack);
    for (int i = 0; i < graph->numVertices; i++) {
        Node* curr = graph->adjlist[i];
        while (curr != NULL) {
            Node* temp = curr;
            curr = curr->next;
            free(temp);
        }
    }
    free(graph->adjlist);
    free(graph->fromTo);
    free(graph);
    
    return count;
}

int main(int argc, char* argv[]) {
    int argv2 = atoi(argv[2]);
    Graph *graph = readGraphFromFile(argv[1], argv2);
    printf("%d\n", dfs(graph));
}
