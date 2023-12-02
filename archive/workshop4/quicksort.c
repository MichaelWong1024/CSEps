#include<stdio.h>
#include<stdlib.h>
#include<string.h>

typedef struct {
    double value1;
    double value2;
    double value3;
} Item;

int quicksortCalls = 0;

void swap(Item* a, Item* b) {
    Item temp = *a;
    *a = *b;
    *b = temp;
}

int partition(Item arr[], int low, int high) {
    double pivot = arr[high].value2;
    int i = low - 1;

    for (int j = low; j <= high - 1; j++) {
        if (arr[j].value2 < pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

void quickSort(Item arr[], int low, int high) {
    quicksortCalls++;
    if (low < high) {
        int pi = partition(arr, low, high);
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

int main(int argc, char *argv[]) {
    FILE *file;
    int N;
    Item* array;

    file = fopen(argv[1], "r");
    if (file == NULL) {
        printf("Error opening file.\n");
        return 1;
    }

    fscanf(file, "%d", &N);
    array = (Item *)malloc(N * sizeof(Item));

    for (int i = 0; i < N; i++) {
        fscanf(file, "%lf %lf %lf", &array[i].value1, &array[i].value2, &array[i].value3);
    }
    fclose(file);

    quickSort(array, 0, N-1);

    if (argc > 2 && strcmp(argv[2], "p") == 0) {
        for (int i = 0; i < N; i++) {
            printf("%.4f\t%.4f\t%.4f\n", array[i].value1, array[i].value2, array[i].value3);
        }
    }

    printf("%d\n", quicksortCalls);

    free(array);
    return 0;
}
