#include <stdio.h>
#include <stdlib.h>

// The struct for linked list nodes.
typedef
struct node {
	int value;
	struct node* next;
}
node;

int main() {
	// The dummy node acts as the initial node for the linked list.
	node* dummy = (node*) malloc(sizeof(node));
	dummy -> value = -1;
	dummy -> next = NULL;

	// A pointer to traverse the list.
	node* currPtr = dummy;

	// Add 10 new nodes with values from 0 to 9 to the front.
	for (int i = 0; i < 10; i++) {
		// Create a new node and allocate memory for it.
		node* newNode = (node*) malloc(sizeof(node));
		newNode -> value = i;

		// Link the new node to the origin node (currPtr).
		newNode -> next = currPtr;
		// Update currPtr to point to the new node.
		currPtr = newNode;
	}

	// Set head to point to the first node in the list.
	node* head = currPtr;
	while (head != NULL) {
		printf("%d", head -> value);
		head = head -> next;
		printf("\n");
	}

	// Reset head to point to the first node again.
	head = currPtr;

	// Free the allocated memory for each node
	while (head != NULL) {
		// Temporarily store the current node
		node* temp = head;

		// Move the pointer to the next node
		head = head -> next;

		// Free the memory of the current temp node
		free(temp);
	}
}
