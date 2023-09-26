// Typedef and function definitions for PQ
#include <stdio.h>
#include <stdlib.h>

// Define a structure for flight information
typedef struct flightinfo {
	char airline; 		// Use U for United, D for Delta, A for American, S for Southwest
	int flightnumber;
	int time;
} Flightinfo;

// Define a structure for items in the priority queue
typedef struct item {
	Flightinfo flight;
	int* priority; // Priority is a pointer to an int, which will point to flight time
	struct item* next;
} Item;

// Initialize a priority queue
Item* Initialize() {
	// Allocate memory for a dummy node to serve as the head of the priority queue
	Item* dummyNode = (Item*)malloc(sizeof(Item));
	if (dummyNode == NULL) {
		return NULL; // Return NULL if memory allocation fails
	}
	// Initialize the dummy node's flight information
	dummyNode->flight.time = -1;
	dummyNode->flight.flightnumber = -1;
	dummyNode->flight.airline = ' ';
	// Set the dummy node's priority to point to the flight time
	dummyNode->priority = &(dummyNode->flight.time);
	dummyNode->next = NULL; // Initialize the next pointer to NULL
	return dummyNode; // Return the initialized priority queue
}

// Add an item to the priority queue
int Add(Item* mypq, char itemAirline, int itemFlightnumber, int itemTime) {
	if (itemFlightnumber <= 0 || itemTime < 0) {
		return 1; // Return an error code if flight number or time is invalid
	}
	Item* newnode = (Item*)malloc(sizeof(Item)); // Allocate memory for a new node
	if (newnode == NULL) {
		return 1; // Return an error code if memory allocation fails
	}
	// Initialize the flight information of the new node
	newnode->flight.time = itemTime;
	newnode->flight.flightnumber = itemFlightnumber;
	newnode->priority = &(newnode->flight.time); // Set priority to point to flight time
	newnode->flight.airline = itemAirline;
	Item* current = mypq; // Initialize a current pointer to the head of the queue
	while (current->next != NULL && *(current->next->priority) <= *(newnode->priority)) {
		// Traverse the queue to find the appropriate position based on priority
		current = current->next;
	}
	newnode->next = current->next; // Insert the new node into the queue
	current->next = newnode;
	return 0; // Return success
}

// Remove the item with the highest priority from the priority queue, if it's not empty
void Remove(Item* mypq) {
	if (mypq == NULL || mypq->next == NULL) {
		return; // Return if the queue is empty or invalid
	}
	Item* temp = mypq->next;
	mypq->next = temp->next;
	free(temp); // Free the memory of the removed item
}

// Print the contents of the priority queue
void Print(Item* mypq) {
	printf("Priority queue contents:\n");
	Item* temp = mypq->next;
	while (temp != NULL) {
		// Print the flight information of each item in the queue
		printf("%c %d %d\n", temp->flight.airline, temp->flight.flightnumber, temp->flight.time);
		temp = temp->next;
	}
}

// Count the number of items in the priority queue
int Count(Item* mypq) {
	int count = 0;
	if (mypq == NULL || mypq->next == NULL) {
		return 0; // Return 0 if the queue is empty or invalid
	} else {
		Item* temp = mypq->next;
		while (temp != NULL) {
			count++; // Increment count for each item in the queue
			temp = temp->next;
		}
		return count; // Return the total count
	}
}

// Count the number of items in the priority queue with a specific airline
int CountAirline(Item* mypq, char myairline) {
	int count = 0;
	if (mypq == NULL || mypq->next == NULL) {
		return 0; // Return 0 if the queue is empty or invalid
	} else {
		Item* temp = mypq->next;
		while (temp != NULL) {
			if (temp->flight.airline == myairline) {
				count++; // Increment count for matching airline
			}
			temp = temp->next;
		}
		return count; // Return the count of items with the specified airline
	}
}

// Count the number of items in the priority queue with flight numbers smaller than a given value
int CountSmaller(Item* mypq, int itemflightnumber) {
	int count = 0;
	if (mypq == NULL || mypq->next == NULL) {
		return 0; // Return 0 if the queue is empty or invalid
	} else {
		Item* temp = mypq->next;
		while (temp != NULL) {
			if (temp->flight.flightnumber < itemflightnumber) {
				count++; // Increment count for smaller flight numbers
			}
			temp = temp->next;
		}
		return count; // Return the count of items with smaller flight numbers
	}
}

// Count the number of items in the priority queue with flight numbers larger than a given value
int CountLarger(Item* mypq, int itemflightnumber) {
	int count = 0;
	if (mypq == NULL || mypq->next == NULL) {
		return 0; // Return 0 if the queue is empty or invalid
	} else {
		Item* temp = mypq->next;
		while (temp != NULL) {
			if (temp->flight.flightnumber > itemflightnumber) {
				count++; // Increment count for larger flight numbers
			}
			temp = temp->next;
		}
		return count; // Return the count of items with larger flight numbers
	}
}

// Count the number of items in the priority queue with priority earlier than a given value
int CountEarlier(Item* mypq, int mypriority) {
	int count = 0;
	if (mypq == NULL || mypq->next == NULL) {
		return 0; // Return 0 if the queue is empty or invalid
	} else {
		Item* temp = mypq->next;
		while (temp != NULL) {
			if (*(temp->priority) < mypriority) {
				count++; // Increment count for earlier priorities
			}
			temp = temp->next;
		}
		return count; // Return the count of items with earlier priorities
	}
}

// Count the number of items in the priority queue with priority later than a given value
int CountLater(Item* mypq, int mypriority) {
	int count = 0;
	if (mypq == NULL || mypq->next == NULL) {
		return 0; // Return 0 if the queue is empty or invalid
	} else {
		Item* temp = mypq->next;
		while (temp != NULL) {
			if (*(temp->priority) > mypriority) {
				count++; // Increment count for later priorities
			}
			temp = temp->next;
		}
		return count; // Return the count of items with later priorities
	}
}

// Finalize the priority queue and free memory
void Finalize(Item* mypq) {
	Item* next;
	Item* curr = mypq->next;
	while (curr != NULL) {
		next = curr->next;
		free(curr);
		curr = next;
	}
	mypq->next = NULL;
}
