#include <stdio.h>
#include <omp.h>

int main() {
    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();

        printf("Thread %d is starting the first task\n", threadID);
        // Simulating a task with a sleep function
        sleep(5); // Each thread sleeps for a different duration

        printf("Thread %d reached the barrier\n", threadID);
        #pragma omp barrier // Threads will wait here until all have reached this point

        printf("Thread %d is starting the second task\n", threadID);
        // Proceed with the second task
    }

    return 0;
}

