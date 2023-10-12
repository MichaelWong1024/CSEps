#include <stdio.h>
#include <stdlib.h>
#include "connections.h"

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("Usage: ./programName <filename> <max_label>\n");
        return 1;
    }

    int maxLabel = atoi(argv[2]);
    if (maxLabel != 1 && maxLabel != 2) {
        printf("Invalid max label. Defaulting to 2.\n");
        maxLabel = 2;
    }

    Graph* graph = readGraphFromFile(argv[1], maxLabel);
    if (!graph) { return 1; }

    printf("%d\n", countComponents(graph, maxLabel));
    freeGraph(graph);
    return 0;
}
