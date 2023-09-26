#include <stdio.h>
#include <stdlib.h>

int get_index(int i, int j, int cols) {
    return i * cols + j;
}

int* load_adjacency_matrix_from_file(int* size, const char* fileName) {
    FILE* file = fopen(fileName, "r");
    if (!file) { 
        perror("Error opening file");
        return NULL; 
    }

    // Reading the matrix size
    if (fscanf(file, "%d", size) != 1) {
        perror("Error reading matrix size");
        fclose(file);
        return NULL;
    }

    int* matrix = malloc((*size) * (*size) * sizeof(int));
    if (!matrix) {
        perror("Memory allocation error");
        fclose(file);
        return NULL;
    }

    // Reading the adjacency matrix
    for (int i = 0; i < *size; i++) {
        for (int j = 0; j < *size; j++) {
            if (fscanf(file, "%d", &matrix[get_index(i, j, *size)]) != 1) {
                perror("Error reading matrix value");
                free(matrix);
                fclose(file);
                return NULL;
            }
        }
    }

    fclose(file);
    return matrix;
}

void display_vertex_degrees(int* matrix, int size) {
    for (int i = 0; i < size; i++) {
        int inDegree = 0;
        int outDegree = 0;
        for (int j = 0; j < size; j++) {
            inDegree += matrix[get_index(j, i, size)];
            outDegree += matrix[get_index(i, j, size)];
        }
        printf("%d: %d, %d\n", i, inDegree, outDegree);
    }
}

void display_two_hop_paths(int *matrix, int size) {
    for (int u = 0; u < size; u++) {
        for (int v = 0; v < size; v++) {
            if (matrix[get_index(u, v, size)]) {
                for (int w = 0; w < size; w++) {
                    if (matrix[get_index(v, w, size)] && w != u) {
                        printf("%d %d %d\n", u, v, w);
                    }
                }
            }
        }
    }
}

int main() {
    const char filename[] = "testgraph.txt";
    int size;
    int* adjacencyMatrix = load_adjacency_matrix_from_file(&size, filename);
    if (adjacencyMatrix) {
        display_vertex_degrees(adjacencyMatrix, size);
        display_two_hop_paths(adjacencyMatrix, size);
        free(adjacencyMatrix);
    } else {
        printf("Error reading the file or allocating memory.\n");
    }
    return 0;
}
