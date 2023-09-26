#include <stdio.h>
#include <math.h>


// The struct that will be returned by method3
typedef
struct singleReturn {
	double value;
	int index;
}
singleReturn;

// A global array storing 100 double values.
double values[100];
	
void method1(double*, int*);
double method2(int*);
singleReturn method3();

int main() {
	// Populate the array.
	for (int i = 0; i < 100; i++) {
		values[i] = 1 + sin(i);
	}

	// Variables to hold max values and corresponding indices.
	double maxValue1;
	int maxValIdx1, maxValIdx2;

	// Results of method1 are stored in maxValue1 and maxValIdx1. (2 pointers)
	method1(&maxValue1, &maxValIdx1);
	printf("Method 1: Max value is %f at index %d\n", maxValue1, maxValIdx1);

    // The result of method2 is stored in maxValue2. (1 pointer, 1 retrun value)
	double maxValue2 = method2(&maxValIdx2);
	printf("Method 2: Max value is %f at index %d\n", maxValue2, maxValIdx2);

    // Results of method3 are stored in a singleReturn struct. (1 struct)
	singleReturn result = method3();
	printf("Method 3: Max value is %f at index %d\n", result.value, result.index);
}

void method1(double* maxValue, int* maxValIdx) {
	double currMaxVal = 0.0;
	int currMaxValIdx = 0;

	// Loop through 'values' to find max value and its index.
	for (int i = 0; i < 100; i++) {
		if (values[i] > currMaxVal) {
			currMaxVal = values[i];
			currMaxValIdx = i;
		}
	}

	// Update the values of the pointers passed as arguments.
	*maxValue = currMaxVal;
	*maxValIdx = currMaxValIdx;
	return;
}

double method2(int* maxValIdx) {
	double currMaxVal = 0.0;
	int currMaxValIdx = 0;
	for (int i = 0; i < 100; i++) {
		if (values[i] > currMaxVal) {
			currMaxVal = values[i];
			currMaxValIdx = i;
		}
	}

	// Update the value of the pointer passed as an argument.
	*maxValIdx = currMaxValIdx;
	return currMaxVal;
}

singleReturn method3() {
	singleReturn single;
	double currMaxVal = 0.0;
	int currMaxValIdx = 0;
	for (int i = 0; i < 100; i++) {
		if (values[i] > currMaxVal) {
			currMaxVal = values[i];
			currMaxValIdx = i;
		}
	}

	// Update the values in the 'singleReturn' struct.
	single.value = currMaxVal;
	single.index = currMaxValIdx;
	return single;
}
