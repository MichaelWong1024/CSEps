#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <number_of_points> <number_of_threads>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
    long long point_count = strtoll(argv[1], NULL, 10);
    int thread_count = atoi(argv[2]);
    
    if (point_count <= 0 || thread_count <= 0) {
        fprintf(stderr, "Error: Inputs must be positive.\n");
        exit(EXIT_FAILURE);
    }
    
    omp_set_num_threads(thread_count);
    
    long long inside = 0;
    double start = omp_get_wtime();
    
    #pragma omp parallel 
    {
        unsigned int seed = (unsigned int)time(NULL) + omp_get_thread_num();
        #pragma omp for reduction(+:inside)
        for (long long i = 0; i < point_count; ++i) {
            double x = rand_r(&seed) / (double)RAND_MAX;
            double y = rand_r(&seed) / (double)RAND_MAX;
            inside += (x * x + y * y <= 1.0);
        }
    }
    
    double pi = 4.0 * inside / point_count;
    double end = omp_get_wtime();
    printf("%f, %f\n", pi, end - start);


    return EXIT_SUCCESS;
}
