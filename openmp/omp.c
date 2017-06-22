/*omp.c
 * C: Joseph Mclaughlin
 * D: Circular Maze Generator (Open MP) for CIS 431, University of Oregon, Spring 2017  
 *  
 */
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define pi 3.141592
#define THETA_WEIGHT 5
#define RADIAL_WEIGHT 40

typedef struct Point {
	float X, Y;
	int Self, Friends[4], EdgeTo[4]; 
	/* Friends[Left, Right, Out, In] */
	/* EdgeTo[i] Contains location in circ_mz of 
	 * edge_t that connects Self <---> Friend[i]
	 */
} node;
typedef struct Edge {
	node *A, *B;
	int Ai, Bi, Weight, Loc, Out;
} edge;

typedef struct CircularMaze {
	int Radius;
	int CurrentRingSize;
	float Delta;
	float * Angles;
	int ** Points;
	int _Size;
	int * Sections;
	int * Connected;
	int NumEdges;
	node ** Nodes; 
	edge * Edges;
	int * Selected;
	int * Queue;
} circular_maze;

node * GetNode(circular_maze * Maze, int Depth) {
	int i = 0;
	while (Depth >= Maze->Sections[i]) {
		Depth -= Maze->Sections[i];
		++i;
	}
	return &Maze->Nodes[i][Depth];
}
int PriorSum(int * Array, int Depth) {
	int i, rv;
	rv = 0;
	if (Depth == 0)
		return rv;
	for (i=0;i<(Depth);++i)
		rv += (Array[i] * 2);
	return rv;
}

typedef struct Link mst_link_t; 
struct Link {
	int V;
	void * Item;
	mst_link_t * Next;
	mst_link_t * Prev;
}; 
typedef struct Queue {
	int Empty, size, *Indices;
	mst_link_t * Tail;
	mst_link_t * Head;
} mst_queue_t;
void Enqueue (mst_queue_t * Queue, mst_link_t * Item) { 
	if (Queue->Empty == 1) {
		Queue->Empty = 0;
		Queue->Head = Item;
		Queue->Tail = Item;
	} 
	else {
		Queue->Tail->Prev = Item;
		Item->Next = Queue->Tail;
		Queue->Tail = Item;
	}
	Queue->size += 1;
}
void IntEnqueue(mst_queue_t * Queue, int Value) {
	mst_link_t * v = malloc(sizeof(mst_link_t));
	v->V = Value;
	Enqueue(Queue, v);
	return;
}
void Dequeue (mst_queue_t * Queue) { 
	if (Queue->size > 0) {	
		mst_link_t * Temp = Queue->Head;
		if (Queue->Head->Prev != NULL) {
			Queue->Head = Queue->Head->Prev;
			free(Temp);
			Queue->size -= 1;

		}
		else if (Queue->Tail != NULL) { 
			Queue->Head = Queue->Tail;
			free(Temp);
			Queue->size -= 1;

		}
		else {
			Queue->Head = NULL;
			Queue->Empty = 1;
			free(Queue->Head);
		}		
	}
}
void PrintLinks (mst_link_t * Link) { 
	printf("%d, ", Link->V);
	if (Link->Prev != NULL)
		PrintLinks(Link->Prev);
	else
		printf("\n");
}	

typedef struct {
	int Length, size;
    mst_link_t * Links;
} mst_heap_t;

void HeapEnqueue (mst_heap_t *h, int Priority, void * Item) {
	int i, j;
    if (h->Length + 1 >= h->size) {
        h->size = h->size ? h->size * 2 : 4;
        h->Links = realloc(h->Links, h->size * sizeof(mst_link_t));
    }
    i = h->Length + 1;
    j = i / 2;
    while (i > 1 && h->Links[j].V > Priority) {
        h->Links[i] = h->Links[j];
        i = j;
        j = j / 2;
    }
    h->Links[i].V = Priority;
    h->Links[i].Item = Item;
    h->Length++;
}
void HeapEnqueue1 (mst_heap_t *h, int Priority, void * Item) {
    int i, Location;
    mst_link_t * Links;

    if (h->Length == 0){
    	h->Length += 1;
    	h->Links[0].V = Priority;
    	h->Links[0].Item = Item;
    	return;
    }

    h->Length += 1;
    #pragma omp critical 
    {
    	Links = malloc(sizeof(mst_link_t) * h->Length * 2 );
    }
    Location = 0;
    Location = h->Length;
    for (i=0;i<h->Length;++i){
    	if (Priority >= h->Links[i].V)
    		Location = i;
    }
    #pragma omp parallel for 
    for (i=0;i<Location;++i){
    	Links[i].V = h->Links[i].V;
    	Links[i].Item = h->Links[i].Item;
    }
    i = Location;
    Links[i].V = Priority;
    Links[i].Item = Item;
    #pragma omp parallel for 
    for (i=Location;i<h->Length;++i){
    	Links[i+1].V = h->Links[i].V;
    	Links[i+1].Item = h->Links[i].Item;
    }

    #pragma omp critical 
    {
    	free(h->Links);
	}
    h->Links = Links;
    return;
}
void HeapEnqueue2 (mst_heap_t *h, int Priority, void * Item) {
    int i, Location;
    if (h->Length + 1 >= h->size) {
        h->size = h->size * 2;
        h->Links = realloc(h->Links, h->size * 2 * sizeof(mst_link_t));
    }

    if (h->Length == 0){ /* initial */
    	h->Length++ ;
    	h->Links[0].V = Priority;
    	h->Links[0].Item = Item;
    	return;
    }
    Location = 0;
    #pragma omp parallel
{
	#pragma omp for schedule(static) /* the important part */
    for (i=0;i<h->Length;++i){ 
    	if (Location == 0 || i < Location) {
	    	if (Priority >= h->Links[i * 2].V){
	    		#pragma omp critical 
	    		{
	    			if (Location==0 || i < Location){
	    				//printf("%d, %d, %d >= %d\n", i, Location, Priority, h->Links[i * 2].V);
	    				Location = i;
	    			}
	    		}
	    	}
	    }
    }

    if (Location == 0)
    	Location = i+1;


    
    #pragma omp for
    for (i=0;i<h->Length;++i){ /* clean up */
    	h->Links[(2 * i) + 1].V = h->Links[2 * i].V;
    	h->Links[(2 * i) + 1].Item = h->Links[2 * i].Item;
    }
    
    #pragma omp for
    for (i=Location+1;i<h->Length+1;++i){ /* shift */
    	h->Links[2 * i].V = h->Links[(2 * i) - 1].V;
    	h->Links[2 * i].Item = h->Links[(2 * i) - 1].Item;
    }
}
	i = Location;
    h->Links[2 * i].V = Priority;
   	h->Links[2 * i].Item = Item;
    h->Length++;

    return;
}
void * HeadDequeue (mst_heap_t *h) {
    int i, j, k;
    void * Item;
    if (!h->Length) {
        return NULL;
    }
    Item = h->Links[1].Item;
    h->Links[1] = h->Links[h->Length];
    h->Length--;
    i = 1;
    while (1) {
        k = i;
        j = 2 * i;
        if (j <= h->Length && h->Links[j].V < h->Links[k].V) {
            k = j;
        }
        if (j + 1 <= h->Length && h->Links[j + 1].V < h->Links[k].V) {
            k = j + 1;
        }
        if (k == i) {
            break;
        }
        h->Links[i] = h->Links[k];
        i = k;
    }
    h->Links[i] = h->Links[h->Length + 1];
    return Item;
}
void * HeadDequeue1 (mst_heap_t *h) {
	h->Length -= 1;
	if (h->Length == -1)
		return NULL;
	return h->Links[h->Length].Item;
}
void * HeadDequeue2 (mst_heap_t *h) {
	h->Length -= 1;
	if (h->Length == -1)
		return NULL;
	return h->Links[h->Length * 2].Item;
}
void * HeadDequeueP (mst_heap_t *h) {
    int i, j, k;
    void * Item;
    if (!h->Length) {
        return NULL;
    }
    Item = h->Links[1].Item;
    h->Links[1] = h->Links[h->Length];
    h->Length--;
    i = 1;
    while (1) {
        k = i;
        j = 2 * i;
        if (j <= h->Length && h->Links[j].V < h->Links[k].V) {
            k = j;
        }
        if (j + 1 <= h->Length && h->Links[j + 1].V < h->Links[k].V) {
            k = j + 1;
        }
        if (k == i) {
            break;
        }
        h->Links[i] = h->Links[k];
        i = k;
    }
    h->Links[i] = h->Links[h->Length + 1];
    return Item;
}
int MazeDequeue(circular_maze * Maze) {
	int i, Weight, rv;
	int * Mins, * Indices;
	int Threads;
	#pragma omp parallel
	{
		Threads = omp_get_num_threads();
	}

	Mins = malloc(sizeof(int) * Threads);
	Indices = malloc(sizeof(int) * Threads);
	for(i=0;i<Threads;++i) {
		Mins[i] = 100;
		Indices[i] = 100;
	} 
	rv = -1;
	Weight = 100;
	#pragma omp parallel for 
	for (i=0;i<(Maze->NumEdges);++i) {
		if (Maze->Queue[i]==1 && Maze->Edges[i].Weight < Mins[omp_get_thread_num()]) {
			Mins[omp_get_thread_num()] = Maze->Edges[i].Weight;
			Indices[omp_get_thread_num()] = i;
		}
	}
	for(i=0;i<Threads;++i) {
		if (Weight > Mins[i]) {
			Weight = Mins[i];
			rv = Indices[i];
		}
	}
	if (Weight == 100 || rv < 0)
		return -1;
	Maze->Queue[rv] = 0;
	#pragma omp critical
	{
		free(Mins);
		free(Indices);
	}
	return rv;
}

void PrintEdge(edge Edge) {
	printf("%f`%f`%f`%f", Edge.A->X, Edge.A->Y, Edge.B->X, Edge.B->Y);
}
void PrintEdges(edge * Edges, int Limit) {
	int i;
	for (i=0;i<Limit;++i)
		printf("((%f,%f),(%f,%f))", Edges[i].A->X, Edges[i].A->Y, Edges[i].B->X, Edges[i].B->Y);
}
void PyPrintEdges(circular_maze * Maze) {
	int i;
	for(i=Maze->NumEdges-1;i>0;--i){
		if(Maze->Selected[i]) {
			PrintEdge(Maze->Edges[i]);
			printf(",");
		}
	}
	PrintEdge(Maze->Edges[0]);
	return;
}
void PrintInfo(double Time, long Iterations, double Dev)
{
	printf("%f:%ld:%f\n", Time, Iterations, Dev);
	return;
}
/*
 * Convert pseudo-polar coords to cartesian coordinates.  
 */
void Cartesify(node * N, int Radius, int ThetaOffest, float ThetaIncrement) {
	float Theta = ThetaOffest * ThetaIncrement * (pi/180.0f);
	N->X = Radius * cos(Theta);
	N->Y = Radius * sin(Theta);
	return;
}
/*
 *
 */
void InitalizeNodes(circular_maze * Maze) {
	int i, j, l, k, Top;
	float Angle;
	for(i=1;i<=Maze->Radius;i++) {
		Angle = (180.0*acos((2.0*(i*i)-1.0)/(2.0*(i*i))))/pi;
		if (Angle*2.0f < Maze->Delta) {
			Maze->CurrentRingSize *= 2;
			Maze->Delta = Angle;
		}
		Maze->Angles[i-1] = Maze->Delta;
		Maze->Points[i-1] = malloc(Maze->CurrentRingSize * sizeof(int)); 
		Top = Maze->_Size + Maze->CurrentRingSize;
		for (j=0;j<Maze->CurrentRingSize;j++)
			Maze->Points[i-1][j] = Maze->_Size + j;
		Maze->_Size = Top;
		Maze->Sections[i-1] = Maze->CurrentRingSize;
	}
	Maze->Connected = malloc(Maze->_Size * sizeof(int));
	for (i=0;i<Maze->_Size;i++)
		Maze->Connected[i] = 0; /* Initialize node represenation to unconnected */
	Maze->Nodes = malloc(Maze->_Size * sizeof(node *));
	k = 0;
	for (i=0;i<Maze->Radius;i++) {
		Maze->Nodes[i] = malloc(Maze->Sections[i] * sizeof(node));
		for (j=0;j<Maze->Sections[i];j++){
			Cartesify(&Maze->Nodes[i][j], i+1, j, (1.0f/Maze->Sections[i])*360.0f);
			Maze->Nodes[i][j].Self = k;
			for (l=0;l<4;l++)
				Maze->Nodes[i][j].Friends[l] = -1;
			++k;
		}
	}
	return;
}
void InitalizeEdges(circular_maze * Maze) {
	int i, j, k, Radius; 
	int n = 1;
	long NumEdges;
	k = 0;
	NumEdges = 0;
	for (i=0;i<Maze->Radius-1;++i)
		NumEdges += Maze->Sections[i] * 2;
	NumEdges += Maze->Sections[i];
	Maze->Edges = malloc(NumEdges * sizeof(edge));
	for (j=0;j<Maze->Sections[0];++j) {
		Maze->Nodes[0][j].Friends[3] = -1;
	}
	Radius = Maze->Radius-1;
	for (i=0;i<Radius;++i) {
		k = PriorSum(Maze->Sections, i);
		for (j=0;j<Maze->Sections[i];++j) {
			/* set "theta" edges. (The ones that make rings.) */
			Maze->Edges[k + j].Weight = ((k-j) % THETA_WEIGHT) + 1;
			Maze->Edges[k + j].Loc = k + j;
			Maze->Edges[k + j].A = &Maze->Nodes[i][j];
			Maze->Edges[k + j].B = &Maze->Nodes[i][(j+1)%Maze->Sections[i]];
			/* set node neighbors */
			Maze->Nodes[i][j].Friends[0] = Maze->Nodes[i][(j+1)%Maze->Sections[i]].Self;
			Maze->Nodes[i][(j+1)%Maze->Sections[i]].Friends[1] = Maze->Nodes[i][j].Self;

			Maze->Nodes[i][j].EdgeTo[0] = k + j; 
			Maze->Nodes[i][(j+1)%Maze->Sections[i]].EdgeTo[1] = k + j;

		}
		if (Maze->Sections[i] < Maze->Sections[i+1])
			n = 2;
		else
			n = 1;
		k += Maze->Sections[i];
		for (j=0;j<Maze->Sections[i];++j) {
			/* set "radial" edges. (The ones that make spokes.)*/
			Maze->Edges[k + j].Weight = ((k+j) % RADIAL_WEIGHT) + 1;
			Maze->Edges[k + j].Loc = k + j;
			Maze->Edges[k + j].A = &Maze->Nodes[i][j];
			Maze->Edges[k + j].B = &Maze->Nodes[i+1][j*n];
			/* set node neighbors */
			Maze->Nodes[i][j].Friends[2] = Maze->Nodes[i+1][j*n].Self;
			Maze->Nodes[i+1][j*n].Friends[3] = Maze->Nodes[i][j].Self;

			Maze->Nodes[i][j].EdgeTo[2] = k + j;
			Maze->Nodes[i+1][j*n].EdgeTo[3] = k + j;
		}
	}
	k = PriorSum(Maze->Sections, i);
	for (j=0;j<Maze->Sections[i];++j) {
		
		Maze->Edges[k + j].Weight = ((k-j) % THETA_WEIGHT) + 1;
		Maze->Edges[k + j].Loc = k + j;
		Maze->Edges[k + j].A = &Maze->Nodes[i][j];
		Maze->Edges[k + j].B = &Maze->Nodes[i][(j+1)%Maze->Sections[i]];
		/* set node neighbors */
		Maze->Nodes[i][j].Friends[0] = Maze->Nodes[i][(j+1)%Maze->Sections[i]].Self;
		Maze->Nodes[i][(j+1)%Maze->Sections[i]].Friends[1] = Maze->Nodes[i][j].Self;
		Maze->Nodes[i][j].Friends[2] = -1; /* Outer ring has no outer friends */

		Maze->Nodes[i][j].EdgeTo[0] = k + j; 
		Maze->Nodes[i][(j+1)%Maze->Sections[i]].EdgeTo[1] = k + j;
	}
	Maze->NumEdges=NumEdges;
	Maze->Queue = malloc(sizeof(int) * NumEdges);
	Maze->Selected = malloc(sizeof(int) * NumEdges);
	for (i=0;i<NumEdges;i++) {
		Maze->Selected[i] = 0;
		Maze->Queue[i] = 0;
	}
	return;
}

void SpanningTree(circular_maze * Maze) {
	int i, Friend, IterArr[5];
	long Total;
	node * Cur;
	edge ** Items;
	double Start, End, Dev;
	mst_heap_t * Queues;
	Total = 0;
	Start = omp_get_wtime(); 
	Items = malloc(4 * sizeof(edge *));
	Queues = malloc(4 * sizeof(mst_heap_t)); 
	for (i=0;i<4;++i) {
		Queues[i].Links = malloc(sizeof(mst_link_t));
		Queues[i].Length = 0;
		Queues[i].size = 0;
	}
	Cur = GetNode(Maze, Maze->_Size/2);
	Maze->Selected[Maze->_Size/2] = 1;
	for (i=0;i<4;++i) {
		Friend = Cur->Friends[i];
		if (Friend != -1 && Maze->Connected[Friend] != 1){
			Maze->Edges[Cur->EdgeTo[i]].Out = Friend;
			HeapEnqueue(&Queues[i], Maze->Edges[Cur->EdgeTo[i]].Weight, &Maze->Edges[Cur->EdgeTo[i]]);
		}
	}
	#pragma omp parallel for schedule(static,1) num_threads(4)
	for (i=0;i<4;++i) {
		int _j, t, FFriend, Iterations;
		node * _Cur;
		Iterations = 0;
		while ((Items[i] = (edge *) HeadDequeue(&Queues[i])) != NULL) {
			Maze->Selected[Items[i]->Loc] = 1;
			_Cur = GetNode(Maze, Items[i]->Out);
			for (_j=0;_j<4;++_j) {
				FFriend = _Cur->Friends[_j];
				if (FFriend > 0) { 
					t = Maze->Connected[FFriend];
					if (t != 1){
						#pragma omp critical
						{
						Maze->Edges[_Cur->EdgeTo[_j]].Out = FFriend;
						Maze->Connected[FFriend] = 1;
						HeapEnqueue(&Queues[i], Maze->Edges[_Cur->EdgeTo[_j]].Weight, &Maze->Edges[_Cur->EdgeTo[_j]]);
						}
					}	
				}
			}
			++Iterations;
		}
		#pragma omp critical
		{	
			IterArr[omp_get_thread_num()] = Iterations;
			Total += Iterations;
		}
		/*printf("%d iterations by thread: %d\n", Iterations, omp_get_thread_num()); */
	}
	#pragma omp critical
	{
		free(Queues[0].Links);
		free(Queues[1].Links);
		free(Queues[2].Links);
		free(Queues[3].Links);
		free(Queues);
		free(Items);
	}
	End = omp_get_wtime();

	for (i=0;i<4;++i) {
		if (IterArr[i] < (Total/4))
			IterArr[i] =  (Total/4) - IterArr[i];
		else
			IterArr[i] =  IterArr[i] - (Total/4);
	}
	Dev = 0.0f;
	for (i=0;i<4;++i)
		Dev += IterArr[i];
	Dev = Dev/4.0f;
	Dev = Dev/(Total/4);
	PrintInfo(End - Start, Total, Dev);
	/* PyPrintEdges(Maze); */
	return;
}
void PrimMST(circular_maze * Maze) {
	int i, Friend, Iterations, Location;
	node * Cur;
	edge * Item;
	double Start, End;
	mst_heap_t * Queue;
	Start = omp_get_wtime();
	Queue = malloc(sizeof(mst_heap_t)); 
	Queue->Links = malloc(sizeof(mst_link_t));
	Queue->Length = 0;
	Queue->size = 0;
	Cur = GetNode(Maze, 0);
	for (i=0;i<4;++i) {
		Friend = Cur->Friends[i];
		if (Friend != -1 && Maze->Connected[Friend] != 1){
			Maze->Edges[Cur->EdgeTo[i]].Out = Friend;
			Maze->Queue[Cur->EdgeTo[i]] = 1;
			/* HeapEnqueue(Queue, Maze->Edges[Cur->EdgeTo[i]].Weight, &Maze->Edges[Cur->EdgeTo[i]]); */
		}
	}
	
	Iterations = 0;
	while ((Location = MazeDequeue(Maze)) != -1) {
		Maze->Selected[Location] = 1;
		Item = &Maze->Edges[Location];
		Cur = GetNode(Maze, Item->Out);
		for (i=0;i<4;++i) {
			Friend = Cur->Friends[i];
			if (Friend > 0) { 
				if (Maze->Connected[Friend] != 1){
					Maze->Edges[Cur->EdgeTo[i]].Out = Friend;
					Maze->Connected[Friend] = 1;
					Maze->Queue[Cur->EdgeTo[i]] = 1;
					/* HeapEnqueue(Queue, Maze->Edges[Cur->EdgeTo[i]].Weight, &Maze->Edges[Cur->EdgeTo[i]]); */
				}
			}
		}
		++Iterations;
	}	
	End = omp_get_wtime();
	printf("%f:%d\n", End-Start, Iterations);
	#pragma omp critical
	{
		free(Queue->Links);
		free(Queue);
	}
	return;
}
void PrimMST1(circular_maze * Maze) {
	int i, Friend, Iterations;
	node * Cur;
	edge * Item;
	double Start, End;
	mst_heap_t * Queue;
	Start = omp_get_wtime();
	Queue = malloc(sizeof(mst_heap_t)); 
	Queue->Links = malloc(sizeof(mst_link_t) * 1);
	Queue->Length = 0;
	Queue->size = Maze->NumEdges;
	Cur = GetNode(Maze, Maze->_Size/2);
	for (i=0;i<4;++i) {
		Friend = Cur->Friends[i];
		if (Friend != -1 && Maze->Connected[Friend] != 1){
			Maze->Edges[Cur->EdgeTo[i]].Out = Friend;
			HeapEnqueue1(Queue, Maze->Edges[Cur->EdgeTo[i]].Weight, &Maze->Edges[Cur->EdgeTo[i]]);
		}
	}
	
	Iterations = 0;
	while (Queue->Length > 0) {
		Item = (edge *) HeadDequeue1(Queue);
		Maze->Selected[Item->Loc] = 1;
		Cur = GetNode(Maze, Item->Out);
		for (i=0;i<4;++i) {
			Friend = Cur->Friends[i];
			if (Friend > 0) { 
				if (Maze->Connected[Friend] != 1){
					Maze->Edges[Cur->EdgeTo[i]].Out = Friend;
					Maze->Connected[Friend] = 1;
					HeapEnqueue1(Queue, Maze->Edges[Cur->EdgeTo[i]].Weight, &Maze->Edges[Cur->EdgeTo[i]]);
				}
			}
		}
		++Iterations;
	}	
	End = omp_get_wtime();
	printf("%f,%d\n", End-Start, Iterations);
	free(Queue->Links);
	free(Queue);
	return;
}
void PrimMST2(circular_maze * Maze) {
	int i, Friend, Iterations;
	node * Cur;
	edge * Item;
	double Start, End;
	mst_heap_t * Queue;
	Start = omp_get_wtime();
	Queue = malloc(sizeof(mst_heap_t)); 
	Queue->Links = malloc(sizeof(mst_link_t) * 2);
	Queue->Length = 0;
	Queue->size = 2;
	Cur = GetNode(Maze, Maze->_Size/2);
	for (i=0;i<4;++i) {
		Friend = Cur->Friends[i];
		if (Friend != -1 && Maze->Connected[Friend] != 1){
			Maze->Edges[Cur->EdgeTo[i]].Out = Friend;
			HeapEnqueue2(Queue, Maze->Edges[Cur->EdgeTo[i]].Weight, &Maze->Edges[Cur->EdgeTo[i]]);
		}
	}
	
	Iterations = 0;
	while ((Item = (edge *) HeadDequeue2(Queue)) != NULL) {
		Maze->Selected[Item->Loc] = 1;
		Cur = GetNode(Maze, Item->Out);
		for (i=0;i<4;++i) {
			Friend = Cur->Friends[i];
			if (Friend > 0 && Maze->Edges[Cur->EdgeTo[i]].Weight > 0) { 
				if (Maze->Connected[Friend] != 1){
					Maze->Edges[Cur->EdgeTo[i]].Out = Friend;
					Maze->Connected[Friend] = 1;
					HeapEnqueue2(Queue, Maze->Edges[Cur->EdgeTo[i]].Weight, &Maze->Edges[Cur->EdgeTo[i]]);
				}
			}
		}
		++Iterations;
	}	
	End = omp_get_wtime();
	printf("%f,%d\n", End-Start, Iterations);
	free(Queue->Links);
	free(Queue);

	//PyPrintEdges(Maze);

	return;
}
int main(int argc, char** argv) {
	int i;
	circular_maze Maze;
    if (argc != 2) {
        printf("Usage: %s <radius>\n",argv[0]);
        return 1;
    }
    Maze.CurrentRingSize 
    				= 6;
    Maze.Radius 	= atoi(argv[1])-1;
    Maze.Delta 		= 360.0f/6.0f;
    Maze.Points 	= malloc(Maze.Radius * sizeof(int *)); 
    Maze.Sections 	= malloc(Maze.Radius * sizeof(int *)); 
    Maze.Angles 	= malloc(Maze.Radius * sizeof(float *)); 
    Maze._Size 		= 0;
    InitalizeNodes(&Maze);
    InitalizeEdges(&Maze);
    PrimMST2(&Maze);
    /*cleanup*/
    #pragma omp critical
	{
	    for (i=0;i<Maze.Radius;i++) {
	    	free(Maze.Points[i]);
	    	free(Maze.Nodes[i]);
	    }
	    free(Maze.Points);
	    free(Maze.Sections);
	    free(Maze.Nodes);
	    free(Maze.Angles);
	    free(Maze.Connected);
	    free(Maze.Edges);
	    free(Maze.Selected);
	    free(Maze.Queue);
	}

    return 0;
}
