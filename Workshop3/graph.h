//  DO NOT modify this file
#include <stdio.h>
#include <stdlib.h>

#define NVERT 7

void getGraph(int *graph) {
    int mygraph[NVERT*NVERT] = {0, 0, 9, 0, 7, 2, 0,
                                0, 0, 5, 4, 0, 0, 6,
                                9, 5, 0, 1, 0, 3, 4,
                                0, 4, 1, 0, 2, 4, 0,
                                7, 0, 0, 2, 0, 0, 1,
                                2, 0, 3, 4, 0, 0, 3,
                                0, 6, 4, 0, 1, 3, 0,};

    for (int i=0; i<NVERT*NVERT; i++) {
        graph[i] = mygraph[i];
    }
    return;
}