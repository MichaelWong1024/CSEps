#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>


// Structure for an EdgeNode: "to" vertex in the adjacency list of a given "from" vertexs
typedef struct EdgeNode{
    int to_vertex; //the "to" vertex number in the adjacency list
    int weight; 
    struct EdgeNode* next; //point to next edgenode 
}EdgeNode;

// Structure for a vertex in the vertex array
typedef struct VertexNode{
    int vertex;
    EdgeNode* firstNeighbor;
}VertexNode;

// Define a structure for a graph
typedef struct{
    VertexNode* adjlist; // an array of vertex
    int n; // Number of vertex nodes
}Graph;

// Function: create a graph
void createGraph(Graph* totalGraph, const char *file_name){
    // Read file
    FILE *fp;
    fp = fopen(file_name, "r");
    if (!fp){
        printf("file read wrong");
        exit(1);
    }

    //Read the number of vertices
    int number;
    fscanf(fp, "%d", &number);
    // printf("%d\n", number);
    totalGraph->n = number;

    // Alocate memory dynamically
    totalGraph->adjlist = (VertexNode*)malloc(sizeof(VertexNode) * number);
    if (totalGraph->adjlist == NULL){
        printf("Error in allocating memory.\n");
        fclose(fp);
        exit(1);
    }

    // Create the vertex node array
    for (int i = 0; i < number; i++){
        totalGraph->adjlist[i].vertex = i;
        totalGraph->adjlist[i].firstNeighbor = NULL;
    }

    // Read edges and store them
    int i, j, classifier;
    EdgeNode* e;

    while(fscanf(fp, "%d %d %d", &i, &j, &classifier) == 3){
        //alocate memory to EdgeNode
        e = (EdgeNode*) malloc(sizeof(EdgeNode));
        if (e == NULL){
            printf("Error in allocating memory.\n");
            fclose(fp);
            exit(1);
        }
        // put vertex j before the adjlist[i]'s current first neighbor and make vertex j become its new neighbor
        e->to_vertex = j;
        e->weight = classifier;
        e->next = totalGraph->adjlist[i].firstNeighbor;
        totalGraph->adjlist[i].firstNeighbor = e;  
    }
    fclose(fp);
}

//Function: use recussion to conduct depthfirstsearch
void DepthFirstSearch(Graph* G, int v, bool visited[], int type){
    //if vertex v is visited, set it to true
    visited[v] = true;

    // traverse all neighbors of v
    EdgeNode* neighbor = G->adjlist[v].firstNeighbor;

    if (type == 2){
        // traverse all contacts(both frequent and occasional)
        while (neighbor != NULL){
            int neighbor_vertex = neighbor->to_vertex;
            // if the neighbor of v is not visited, traverse it
            if (!visited[neighbor_vertex]){
                // printf("The vertex number you have just visited is: %d\n", neighbor_vertex);
                DepthFirstSearch(G, neighbor_vertex, visited, type);
            }
            neighbor = neighbor->next;
        }
    }else{
        // traverse frequent contacts only
        while (neighbor != NULL){
            int neighbor_vertex = neighbor -> to_vertex;
            int neighbor_weight = neighbor -> weight;
            // if the neighbor of v is not visited and the weight is 1, traverse it
            if (!visited[neighbor_vertex] && neighbor_weight == 1){
                // printf("The vertex number you have just visited is: %d\n", neighbor_vertex);
                DepthFirstSearch(G, neighbor_vertex, visited, type);
            }
            neighbor = neighbor->next;
        }
    }
}

void free_func(Graph *graph){
    if (graph == NULL){
        return;
    }
    for (int i = 0; i<graph->n; i++){
        EdgeNode* temp = graph->adjlist[i].firstNeighbor;
        while (temp != NULL){
            EdgeNode* next = temp->next;
            free(temp);
            temp = next;
        }
    }
    free(graph->adjlist);
    graph->n = 0;
    graph->adjlist = NULL;
}

int main(int argc, char *argv[]){
    //check if the input file name is correct
    const char* file_name;
    if(argc >=2){
        file_name = argv[1];
    }else{
        printf("Wrong file name");
        return 1;
    }
    
    Graph totalGraph;
    createGraph(&totalGraph, file_name);

    int numberOfVertex = totalGraph.n;
    bool visited[numberOfVertex];

    for (int i = 0; i < numberOfVertex; i++){
        visited[i] = false;
    }

    // Conduct DFS
    int count = 0;
    // Check what kind of social network we should search
    // argument[3]: 1, only connections with label of 1; 2, connections with labels of 1 or 2; otherwise, choose the closest one
    int type;
    if(strcmp(argv[2], "1") == 0){
        type = 1;
    }else if(strcmp(argv[2], "2") == 0){
        type = 2;
    }else{
        double input_number = atof(argv[2]);
        double diff_1 = abs(input_number-1.0);
        double diff_2 = abs(input_number-2.0);
        if (diff_1 > diff_2){
            type = 2;
        }else{
            type = 1;
        }
    } 

    for(int i = 0; i < numberOfVertex; i++){
        if (!visited[i]){
            // printf("The unvisited vertex number is: %d\n", i);
            DepthFirstSearch(&totalGraph, i, visited, type);
            count++;
        }
    }
    printf("%d\n", count);

    // free memory
    free_func(&totalGraph);
}