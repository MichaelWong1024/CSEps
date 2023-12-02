#include <stdio.h>
#include <math.h>

#define MAX_ITERATIONS 50 // // Maximum number of iterations
#define TOLERANCE 1e-10 // Convergence tolerance

/*
status:
    0: success
    1: max iterations reached
    2: singular jacobian
    3: negative in sqrt
    4: not implemented
*/
enum Status {
    SUCCESS,
    MAX_ITERATIONS_REACHED,
    SINGULAR_JACOBIAN,
    NEGATIVE_IN_SQRT,
    NOT_IMPLEMENTED
};

double f(double x, double y);
double g(double x, double y);
void jacobian_calculation(double x, double y, double j[2][2], double* det);
int has_converged(double x, double y, double newX, double newY);


void newton_method(double x0, double y0, double* x, double* y, int* iterations, int* status) {
    *x = x0; // initial guess for x
    *y = y0; // initial guess for y
    *iterations = 0;
    *status = NOT_IMPLEMENTED;

    double j[2][2], det; // Jacobian matrix and its determinant
    double newX, newY; // New values of x and y after each iteration

    while (*iterations < MAX_ITERATIONS) {
        double currX = *x, currY = *y;

        jacobian_calculation(currX, currY, j, &det); // Calculate Jacobian matrix and its determinant

        // Check for a singular Jacobian matrix
        if (det == 0) {
            *status = SINGULAR_JACOBIAN;
            return;
        }

        // Update x and y using the Newton's method
        newX = currX -(j[1][1] * f(currX, currY) - j[0][1] * g(currX, currY)) / det;
        newY = currY -(-j[1][0] * f(currX, currY) + j[0][0] * g(currX, currY)) / det;

        // Check for convergence
        if (has_converged(currX, currY, newX, newY)) {
            *status = SUCCESS;
            return;
        }

        // Update x and y for the next iteration
        *x = newX;
        *y = newY;
        (*iterations)++;

    }
    *status = MAX_ITERATIONS_REACHED;
    return;
}

void fixed_point_iteration(double x0, double y0, double* x, double* y, int* iterations, int* status) {
    *x = x0;
    *y = y0;
    *iterations = 0;
    *status = NOT_IMPLEMENTED;

    while (*iterations < MAX_ITERATIONS) {
	    double currX = *x, currY = *y;

        // Check for negative values before taking square roots
	    if (currX < 0 || currY < 0) {
		    *status = NEGATIVE_IN_SQRT;
		    return;
	    }

        // Update x and y using the Fixed-point iteration
	    double newX = sqrt(1 - currY * currY), newY = sqrt(9 - 5 * currX * currX) / sqrt(21);

        // Check for convergence
	    if (has_converged(currX, currY, newX, newY)) {
            	*status = SUCCESS;
           	 return;
	    }

        // Update x and y for the next iteration
	    *x = newX;
	    *y = newY;
	    (*iterations)++;
    }
    *status = MAX_ITERATIONS_REACHED;
    return; 
}

double f(double x, double y) {
    return x * x + y * y - 1;
}

double g(double x, double y) {
    return 5 * x * x + 21 * y * y - 9;
}

void jacobian_calculation(double x, double y, double j[2][2], double* det) {
	j[0][0] = 2 * x;
	j[0][1] = 2 * y;
	j[1][0] = 10 * x;
	j[1][1] = 42 * y;
	*det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
}

int has_converged(double currX, double currY, double newX, double newY) {
	return fabs(currX - newX) < TOLERANCE && fabs(currY - newY) < TOLERANCE;
}
