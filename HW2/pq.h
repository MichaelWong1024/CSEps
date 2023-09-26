// Item and function prototypes for priority queue
typedef struct item Item;

Item* Initialize();
int Add(Item* mypq, char itemAirline, int itemFlightnumber, int itemTime);
void Remove(Item* mypq);
void Print(Item* mypq);
int Count(Item* mypq);
int CountAirline(Item* mypq, char myairline);
int CountSmaller(Item* mypq, int itemflightnumber);
int CountLarger(Item* mypq, int itemflightnumber);
int CountEarlier(Item* mypq, int mypriority);
int CountLater(Item* mypq, int mypriority);
void Finalize(Item* mypq);
