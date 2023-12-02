#include <stdio.h>
#include <omp.h>

int main() {
    // Starting a parallel region
    #pragma omp parallel
    {
        int threadID = omp_get_thread_num(); // Get the ID of the current thread
        printf("Hello from thread %d\n", threadID); // Each thread prints its ID
    }
    // End of the parallel region

    return 0;
}

