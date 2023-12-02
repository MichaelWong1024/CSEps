// Priority queue using a linked list
#include <stdio.h>
#include <stdlib.h>
#include "pq.h"

int main()
{
    Item* mypq = Initialize();

    // Test cases with an empty queue
    printf("Initialization finished\n");
    Print(mypq);
    printf("Number of flights: %d\n", Count(mypq));
    Remove(mypq);
    Remove(mypq);
    printf("Remove process finished\n");
    Print(mypq);
    printf("Number of flights: %d\n\n", Count(mypq));

    printf("Start Adding:\n");

    // Test cases without errors
    int returnval;
    returnval = Add(mypq, 'S', 1402, 2317);
    returnval = Add(mypq, 'V', 0, 0);
    if (returnval == 1) printf("Error adding to pq\n");
    returnval = Add(mypq, 'C', 70, 1043);
    returnval = Add(mypq, 'S', 201, 500);
    returnval = Add(mypq, 'D', 1691, 1601);
    returnval = Add(mypq, 'X', 1500, 1300);
    returnval = Add(mypq, 'Y', 10, 130);

    Print(mypq);

    printf("Number of flights: %d\n", Count(mypq));
    printf("Number of Southwest flights: %d\n", CountAirline(mypq, 'S'));
    printf("Number of Delta flights: %d\n", CountAirline(mypq, 'D'));

    int flighttest = 100;
    printf("Number of flight numbers smaller/larger than %d: %d, %d\n", flighttest,
           CountSmaller(mypq, flighttest), CountLarger(mypq, flighttest));

    int mypriority = 1200;
    printf("Number of flights earlier/later than %d: %d, %d\n", mypriority,
           CountEarlier(mypq, mypriority), CountLater(mypq, mypriority));
    Remove(mypq);
    Remove(mypq);
    printf("Remove Processed finished\n");
    Print(mypq);

    // Additional Test Cases
    printf("\n--- Additional Test Cases ---\n");

    // Adding duplicate flight numbers
    Add(mypq, 'A', 123, 1000);
    Add(mypq, 'B', 123, 1000);
    Print(mypq);

    // Large number of flights
    for (int i = 0; i < 1000; i++) {
        Add(mypq, 'A', i, 1000 + i);
    }
    printf("After adding 1000 flights: Number of flights: %d\n", Count(mypq));

    // Handling invalid airlines
    Add(mypq, '?', 1234, 1200);
    Print(mypq);

    // Priority ordering
    Add(mypq, 'A', 200, 700);
    Add(mypq, 'A', 201, 800);
    Add(mypq, 'A', 202, 600);
    Print(mypq);

    // Add and remove continuously
    Add(mypq, 'A', 250, 1300);
    Remove(mypq);
    Add(mypq, 'A', 251, 1400);
    Remove(mypq);
    Print(mypq);

    // Query on empty list
    Finalize(mypq);
    printf("Number of Delta flights after emptying: %d\n", CountAirline(mypq, 'D'));

    // Add flights with decreasing time
    for (int i = 0; i < 10; i++) {
        Add(mypq, 'A', i, 1000 - i * 10);
    }
    Print(mypq);

    // Add flights with the same time but different airlines
    Add(mypq, 'A', 300, 1300);
    Add(mypq, 'B', 301, 1300);
    Print(mypq);

    Finalize(mypq);
    free(mypq);

    return 0;
}
