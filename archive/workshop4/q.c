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