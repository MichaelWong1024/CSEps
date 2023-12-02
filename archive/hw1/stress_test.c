#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "solution.h"

// Tolerance level for floating-point comparisons
#define TOLERANCE 1e-10

enum Status {
    SUCCESS,
    MAX_ITERATIONS_REACHED,
    SINGULAR_JACOBIAN,
    NEGATIVE_IN_SQRT,
    NOT_IMPLEMENTED
};

// Expected intersection points (x, y)
double answers[4][2] = {
    {sqrt(3)/2.0, 1/2.0},
    {sqrt(3)/2.0, -1/2.0},
    {-sqrt(3)/2.0, 1/2.0},
    {-sqrt(3)/2.0, -1/2.0}
};

void test_newton_method(double x, double y);
void test_fixed_point_iteration(double x, double y);

int main() {
    for (double x = -10000; x < 10000; x += 111.1) {
        for (double y = -10000; y < 10000; y += 111.1) {
            test_newton_method(x, y);
        }
    }

    for (double x = -100; x < 100; x += 11.11) {
	    for (double y = -100; y < 100; y += 11.11) {
		    test_fixed_point_iteration(x, y);
	    }
    }
    
    printf("All tests passed.\n");
    return 0;
}

void test_newton_method(double x0, double y0) {
    double x, y;
    int iterations, status;
    newton_method(x0, y0, &x, &y, &iterations, &status);

    // If the method didn't converge, exit this test iteration
    switch(status) {
        case SINGULAR_JACOBIAN:
            printf("Method failed because of the singular jacobian matrix\n");
            return;
        case MAX_ITERATIONS_REACHED:
            printf("Max iterations reached\n");
            return;
        case SUCCESS:
            break;  // Continue on to test for known intersection points
        default:
            return;  // Not implemented or other unknown status
    }

    // Loop through known intersection points to see if calculated (x, y) matches
    for (int i = 0; i < 4; i++) {
        if (fabs(x - answers[i][0]) < TOLERANCE && fabs(y - answers[i][1]) < TOLERANCE) {
            return;
        }
    }
    // If none of the known points match, assert failure
    assert(0);
}

void test_fixed_point_iteration(double x0, double y0) {
    double x, y;
    int iterations, status;
    fixed_point_iteration(x0, y0, &x, &y, &iterations, &status);

    switch(status) {
        case NEGATIVE_IN_SQRT:
            printf("x = %f, y = %f, one of them is negative\n", x, y);
            return;
        case MAX_ITERATIONS_REACHED:
            printf("Max iterations reached\n");
            return;
        case SUCCESS:
            break;  // Continue on to test for known intersection points
        default:
            return;  // Not implemented or other unknown status
    }

    // Loop through known intersection points to see if calculated (x, y) matches
    for (int i = 0; i < 4; i++) {
        if (fabs(x - answers[i][0]) < TOLERANCE && fabs(y - answers[i][1]) < TOLERANCE) {
            return;
        }
    }
    // If none of the known points match, assert failure
    assert(0);
}
