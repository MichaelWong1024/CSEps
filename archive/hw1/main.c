#include <stdio.h>

#include "solution.h"

int main() {
    double x0, y0, x_newton, y_newton, x_fixed_point, y_fixed_point;
    int iterations_newton, status_newton, iterations_fixed_point, status_fixed_point;

    // take input using keyboard
    printf("Enter x0: ");
    scanf("%lf", &x0);
    printf("Enter y0: ");
    scanf("%lf", &y0);
    printf("\n");

    // call newton method
    newton_method(x0, y0, &x_newton, &y_newton, &iterations_newton, &status_newton);

    // print results
    printf("Newton method:\n");
    printf("x: %.10lf, y: %.10lf, iterations: %d, status: %d\n", x_newton, y_newton, iterations_newton, status_newton);

    // call fixed point iteration
    fixed_point_iteration(x0, y0, &x_fixed_point, &y_fixed_point, &iterations_fixed_point, &status_fixed_point);

    // print results
    printf("Fixed point iteration:\n");
    printf("x: %.10lf, y: %.10lf, iterations: %d, status: %d\n", x_fixed_point, y_fixed_point, iterations_fixed_point, status_fixed_point);
    return 0;
}