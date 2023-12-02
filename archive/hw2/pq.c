// Typedef and function definitions for PQ
#include <stdio.h>
#include <stdlib.h>

// Define Flightinfo struct to hold flight information
typedef struct flightinfo {
	char airline; 		// use U for United, D for Delta, A for American, S for Southwest
	int flightnumber;
	int time; 		// Departure time in military format
} Flightinfo;

// Define Item struct to represent a Priority Queue node
typedef struct item {
	Flightinfo flight; 		// The flight information
	int* priority; // priority is a pointer to an int, will point to flighttime
	struct item * next;
} Item;

// Initializes Priority Queue with a dummy node
// Returns pointer to the dummy node
Item* Initialize()
{
	Item* dummy = malloc(sizeof(Item));
	if (dummy == NULL) { return NULL; }
	Flightinfo dummyFlight = {
    .airline = ' ',
    .flightnumber = -1,
    .time = -1
	};
	dummy->flight = dummyFlight;
	dummy->priority = &(dummy->flight.time);
	dummy->next = NULL;
	return dummy;
}

// Adds a flight to the Priority Queue
// Returns 0 if successful, 1 otherwise
int Add(Item* mypq, char itemAirline, int itemFlightnumber, int itemTime)
{
	if (itemFlightnumber <= 0 || itemTime < 0) { return 1; }
	Item* newItem = malloc(sizeof(Item));
	Flightinfo newFlight = {
    .airline = itemAirline,
    .flightnumber = itemFlightnumber,
    .time = itemTime
	};
	newItem->flight = newFlight;
	newItem->priority = &(newItem->flight.time);

	Item* current = mypq;
	while (current->next != NULL && *(newItem->priority) > *(current->next->priority)) {
		current = current->next;
	}
	newItem->next = current->next;
	current->next = newItem;
	return 0;
}

// Removes the highest priority flight from the Priority Queue
void Remove(Item* mypq)
{
	if (mypq->next == NULL) { return; }

	Item* headToBeFreed = mypq->next;
	mypq->next = mypq->next->next;
	free(headToBeFreed);
}

// Prints the flights in the Priority Queue
void Print(Item* mypq)
{
	printf("Priority queue contents:\n");

	Item* current = mypq->next;
	while (current != NULL) {
		printf("%c %d %d\n", current->flight.airline, current->flight.flightnumber, current->flight.time);
		current = current->next;
	}
}

// Counts the number of flights in the Priority Queue
// Returns the count
int Count(Item* mypq)
{
	int count = 0;
	Item* current = mypq->next;
	while (current != NULL) {
		count++;
		current = current->next;
	}
	return count;
}

// Counts the number of flights for a specific airline
// Returns the count
int CountAirline(Item* mypq, char myairline)
{
	int count = 0;
	Item* current = mypq->next;
	while (current != NULL) {
		if (current->flight.airline == myairline) { count++; }
		current = current->next;
	}
	return count;
}

// Counts the number of flights with flight numbers smaller than a given flight number
// Returns the count
int CountSmaller(Item* mypq, int itemflightnumber)
{
	int count = 0;
	Item* current = mypq->next;
	while (current != NULL) {
		if (current->flight.flightnumber < itemflightnumber) { count++; }
		current = current->next;
	}
	return count;
}

// Counts the number of flights with flight numbers larger than a given flight number
// Returns the count
int CountLarger(Item* mypq, int itemflightnumber)
{
	int count = 0;
	Item* current = mypq->next;
	while (current != NULL) {
		if (current->flight.flightnumber > itemflightnumber) { count++; }
		current = current->next;
	}
	return count;
}

// Counts the number of flights departing earlier than a given time
// Returns the count
int CountEarlier(Item* mypq, int mypriority)
{
	int count = 0;
	Item* current = mypq->next;
	while (current != NULL) {
		if (*(current->priority) < mypriority) { count++; }
		current = current->next;
	}
	return count;
}

// Counts the number of flights departing later than a given time
// Returns the count
int CountLater(Item* mypq, int mypriority)
{
	int count = 0;
	Item* current = mypq->next;
	while (current != NULL) {
		if (*(current->priority) > mypriority) { count++; }
		current = current->next;
	}
	return count;
}

// Frees all the nodes in the Priority Queue
void Finalize(Item* mypq)
{
    Item* temp = mypq->next;
    Item* toBeRemoved;
    while (temp != NULL) {
        toBeRemoved = temp;
        temp = temp->next;
        free(toBeRemoved);
    }
    // Reset the next pointer of dummy node to NULL
    mypq->next = NULL;
}

