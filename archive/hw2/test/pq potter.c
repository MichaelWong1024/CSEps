// Typedef and function definitions for PQ
#include <stdio.h>
#include <stdlib.h>

typedef struct flightinfo {
	char airline; 		// use U for United, D for Delta, A for American, S for Southwest
	int flightnumber;
	int time;
} Flightinfo;

typedef struct item {
	Flightinfo flight;
	int* priority; // priority is a pointer to an int, will point to flighttime
	struct item* next;
} Item;

/// @brief Initialize a dummy variable that will be the head of the priority queue;
/// @return return the dummy variable.
Item* Initialize()
{ 
    Item* dummy = (Item*)malloc(sizeof(Item));
    if (dummy == NULL) {
        // Handle memory allocation error
        return NULL;
    }
	dummy -> flight.airline =' ';
	dummy -> flight.flightnumber = -1;
	dummy -> flight.time =-1;
	dummy -> priority = &(dummy->flight.time);
	dummy -> next = NULL;
	return dummy;
}

/// @brief it first check whether a non-positive value was entered for the itemFlightnumber or a 
/// negative value for the itemTime. Make a node for the filigh and add it into the queue
/// @return Return 0 if the operation succeeds, or 1 if it fails
int Add(Item* mypq, char itemAirline, int itemFlightnumber, int itemTime)
{
	if(itemFlightnumber<=0 ||itemTime <0){
		return 1;
	}
	Item* newItem = (Item*)malloc(sizeof(Item));
	Item* curr = mypq;
	Item* temp;
    if (newItem == NULL) {
        // Handle memory allocation error
        return 1;
    }
	newItem -> flight.airline =itemAirline;
	newItem -> flight.flightnumber = itemFlightnumber;
	newItem -> flight.time =itemTime;
	newItem -> priority = &(newItem->flight.time);
	newItem -> next = NULL;
	/// @brief curr starts from the begining. its next is not NULL, it starts to add;
	/// if next node's priority is bigger, we add the node 
	/// else, curr move forward
	while(curr->next!=NULL){
		if(*(curr->next->priority)>*(newItem->priority)){
			temp = curr->next;
			curr->next  = newItem;
			newItem->next = temp;
			return 0;
		}
		curr = curr->next;
	}
	curr->next=newItem;
	return 0;
}

void Remove(Item* mypq){
	Item* temp;
	if(mypq->next != NULL){
		temp = mypq->next;
		mypq->next = temp->next;
		free(temp);
	}
}


void Print(Item* mypq)
{
	Item* temp;
	printf("Priority queue contents:\n");
	temp = mypq;
	/// @brief temp points to the head node and keep moving forward.
	while(temp != NULL){
		if(*(temp->priority)!= -1){
			printf("%c %d %d\n",temp->flight.airline,
			temp->flight.flightnumber, temp->flight.time);
		}
		temp = temp->next;
	}
	

}
/// @brief Count the number of airlines in the queue
int Count(Item* mypq)
{
	int count = 0;
	Item* temp;
	temp =mypq->next;
	while(temp != NULL){
		count=count + 1;
		temp = temp->next;
	}
	return count;
}
/// @brief Count the number of airlines match the name in the queue
int CountAirline(Item* mypq, char myairline)
{
	int count = 0;
	Item* temp;
	temp =mypq->next;
	while(temp != NULL){
		if(temp->flight.airline == myairline){
			count=count + 1;
		}
		temp = temp->next;
	}
	return count;
}
/// @brief Count the number of airlines that is smaller than the provided flight number in the queue
int CountSmaller(Item* mypq, int itemflightnumber)
{
	int count = 0;
	Item* temp;
	temp =mypq->next;
	while(temp != NULL){
		if(temp->flight.flightnumber <itemflightnumber){
			count=count + 1;
		}
		temp = temp->next;
	}
	return count;
}
/// @brief Count the number of airlines that is larger than the provided flight number in the queue
int CountLarger(Item* mypq, int itemflightnumber)
{
	int count = 0;
	Item* temp;
	temp =mypq->next;
	while(temp != NULL){
		if(temp->flight.flightnumber >itemflightnumber){
			count=count + 1;
		}
		temp = temp->next;
	}
	return count;
}
/// @brief Count the number of airlines that is earlier than the provided flight time in the queue
int CountEarlier(Item* mypq, int mypriority)
{
	int count = 0;
	Item* temp;
	temp =mypq->next;
	while(temp != NULL){
		if(*(temp->priority) <mypriority){
			count=count + 1;
		}
		temp = temp->next;
	}
	return count;
}
/// @brief Count the number of airlines that is later than the provided flight time in the queue
int CountLater(Item* mypq, int mypriority)
{
	int count = 0;
	Item* temp;
	temp =mypq->next;
	while(temp != NULL){
		if(*(temp->priority) >mypriority){
			count=count + 1;
		}
		temp = temp->next;
	}
	return count;
}
/// @brief it frees all the memory occupied.
void Finalize(Item* mypq)
{
	Item* next;
	Item* curr = mypq->next;
	while(curr!=NULL){
		next = curr->next;
		free(curr);
		curr = next;	
	}
	mypq->next = NULL;
}


