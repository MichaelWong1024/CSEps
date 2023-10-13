#include <stdio.h>
#include <stdlib.h>
#include "connected.h"

int main(int argc, char* argv[]) {
    // Check if the correct number of command line arguments is provided.
    if (argc != 3) {
        printf("Usage: ./programName <filename> <max_label>\n");
        return 1; // Exit with an error code.
    }

    // Parse the max label from command line arguments.
    int maxLabel = atoi(argv[2]);

    // Validate the provided max label value.
    if (maxLabel != 1 && maxLabel != 2) {
        printf("Invalid max label. Defaulting to 2.\n");
        maxLabel = 2; // Default to 2 if an invalid value is provided.
    }

    // Read the graph from the provided file.
    Graph* graph = readGraphFromFile(argv[1], maxLabel);
    if (!graph) { 
        return 1; // Exit with an error code if reading fails.
    }

    printf("%d\n", countComponents(graph, maxLabel));

    // Deallocate memory used by the graph.
    freeGraph(graph);

    return 0;
}
