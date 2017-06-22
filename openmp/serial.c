/*l.c
 * C: Joseph Mclaughlin
 * D: Circular Maze Generator (serial) for CIS 431, University of Oregon, Spring 2017  
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

typedef struct Link mst_link_t; 
struct Link {
	int V;
	void * Item;
	mst_link_t * Next;
	mst_link_t * Prev;
}; 
typedef struct Queue {
	int Empty, size;
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
		else if (Queue->Tail != NULL) { /*TODO: NOT SAFE AT ALL!*/
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
    printf("Added %d\n", Priority);
    for (i=0;i<h->size;++i)
    	printf("%d, ", h->Links[i].V);
    printf("\n\n");
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
    Links = malloc(sizeof(mst_link_t) * h->Length * 2);
    Location = 0;
    Location = h->Length;
    for (i=0;i<h->Length;++i){
    	if (Priority >= h->Links[i].V)
    		Location = i;
    }
    for (i=0;i<Location;++i){
    	Links[i].V = h->Links[i].V;
    	Links[i].Item = h->Links[i].Item;
    }
    i = Location;
    Links[i].V = Priority;
    Links[i].Item = Item;
    for (i=Location;i<h->Length;++i){
    	Links[i+1].V = h->Links[i].V;
    	Links[i+1].Item = h->Links[i].Item;
    }

    
    free(h->Links);
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
    for (i=0;i<h->Length;++i){ 
    	if (Priority >= h->Links[i * 2].V){
    		Location = i;
    		break;
    	}
    	Location = i+1;
    }	
    for (i=0;i<h->Length;++i){ /* clean up */
    	h->Links[(2 * i) + 1].V = h->Links[2 * i].V;
    	h->Links[(2 * i) + 1].Item = h->Links[2 * i].Item;
    }
    for (i=Location+1;i<h->Length+1;++i){
    	h->Links[2 * i].V = h->Links[(2 * i) - 1].V;
    	h->Links[2 * i].Item = h->Links[(2 * i) - 1].Item;
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

int MazeDequeue(circular_maze * Maze) {
	int i, Weight, rv;
	rv = -1;
	Weight = 100;
	for (i=0;i<(Maze->NumEdges);++i) {
		if (Maze->Queue[i]==1 && Maze->Edges[i].Weight < Weight) {
			Weight = Maze->Edges[i].Weight;
			rv = i;
		}
	}
	if (Weight == 100 || rv < 0)
		return -1;
	Maze->Queue[rv] = 0;
	return rv;
}

void PrintEdge(edge Edge) {
	printf("%f`%f`%f`%f", Edge.A->X, Edge.A->Y, Edge.B->X, Edge.B->Y);
}
void PrintEdges(edge * Edges, int Limit) {
	int i;
	for (i=0;i<Limit;++i)
		printf("((%f,%f),(%f,%f)),\\\n", Edges[i].A->X, Edges[i].A->Y, Edges[i].B->X, Edges[i].B->Y);
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
/*
 * Convert pseudo-polar coords to cartesian coordinates.  
 */
void Cartesify(node * N, int Radius, int ThetaOffest, float ThetaIncrement) {
	float Theta = ThetaOffest * ThetaIncrement * (pi/180.0f);
	N->X = Radius * cos(Theta);
	N->Y = Radius * sin(Theta);
	return;
}

void InitalizeNodes(circular_maze * Maze) {
	int i, j, Top, l, k;
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
	int i, j, k; 
	int n = 1;
	long NumEdges = 0;
	NumEdges = 0;
	for (i=0;i<Maze->Radius-1;++i)
		NumEdges += Maze->Sections[i] * 2;
	NumEdges += Maze->Sections[i];
	Maze->Edges = malloc(NumEdges * sizeof(edge));
	k = 0;
	for (j=0;j<Maze->Sections[0];++j)
		Maze->Nodes[0][j].Friends[3] = -1;
	for (i=0;i<(Maze->Radius-1);++i) {
		for (j=0;j<Maze->Sections[i];++j) {
			/* set "theta" edges. (The ones that make rings.) */
			Maze->Edges[k].Weight = ((k-j) % THETA_WEIGHT) + 1;
			Maze->Edges[k].Loc = k;
			Maze->Edges[k].A = &Maze->Nodes[i][j];
			Maze->Edges[k].B = &Maze->Nodes[i][(j+1)%Maze->Sections[i]];
			/* set node neighbors */
			Maze->Nodes[i][j].Friends[0] = Maze->Nodes[i][(j+1)%Maze->Sections[i]].Self;
			Maze->Nodes[i][(j+1)%Maze->Sections[i]].Friends[1] = Maze->Nodes[i][j].Self;
			Maze->Nodes[i][j].EdgeTo[0] = k; 
			Maze->Nodes[i][(j+1)%Maze->Sections[i]].EdgeTo[1] = k;
			++k;
		}
		if (Maze->Sections[i] < Maze->Sections[i+1])
			n = 2;
		else
			n = 1;
		for (j=0;j<Maze->Sections[i];++j) {
			/* set "radial" edges. (The ones that make spokes.)*/
			Maze->Edges[k].Weight = ((k+j) % RADIAL_WEIGHT) + 1;
			Maze->Edges[k].Loc = k;
			Maze->Edges[k].A = &Maze->Nodes[i][j];
			Maze->Edges[k].B = &Maze->Nodes[i+1][j*n];
			/* set node neighbors */
			Maze->Nodes[i][j].Friends[2] = Maze->Nodes[i+1][j*n].Self;
			Maze->Nodes[i+1][j*n].Friends[3] = Maze->Nodes[i][j].Self;
			Maze->Nodes[i][j].EdgeTo[2] = k;
			Maze->Nodes[i+1][j*n].EdgeTo[3] = k;
			++k;
		}
	}
	for (j=0;j<Maze->Sections[i];++j) {
		Maze->Edges[k].Weight = ((k-j) % THETA_WEIGHT) + 1;
		Maze->Edges[k].Loc = k;
		Maze->Edges[k].A = &Maze->Nodes[i][j];
		Maze->Edges[k].B = &Maze->Nodes[i][(j+1)%Maze->Sections[i]];
		/* set node neighbors */
		Maze->Nodes[i][j].Friends[0] = Maze->Nodes[i][(j+1)%Maze->Sections[i]].Self;
		Maze->Nodes[i][(j+1)%Maze->Sections[i]].Friends[1] = Maze->Nodes[i][j].Self;
		Maze->Nodes[i][j].Friends[2] = -1; /* Outer ring has no outer friends */
		Maze->Nodes[i][j].EdgeTo[0] = k; 
		Maze->Nodes[i][(j+1)%Maze->Sections[i]].EdgeTo[1] = k;
		++k;
	}
	Maze->NumEdges=NumEdges;
	Maze->Selected = malloc(sizeof(int) * Maze->NumEdges);
	Maze->Queue = malloc(sizeof(int) * Maze->NumEdges);
	for (i=0;i<Maze->NumEdges;++i) {
		Maze->Selected[i] = 0;
		Maze->Queue[i] = 0;
	}
	return;
}
void SpanningTree(circular_maze * Maze) {
	int i, Friend;
	node * Cur;
	mst_queue_t * Queue = malloc(sizeof(mst_queue_t)); 
	Queue->Empty = 1; 
	Queue->size = 0;
	Cur = GetNode(Maze, Maze->_Size/2);
	for (i=0;i<4;++i) {
		Friend = Cur->Friends[i];
		if (Friend != -1 && Maze->Connected[Friend] != 1){
			mst_link_t * l = malloc(sizeof(mst_link_t));
			l->V = Friend;
			l->Item = &Maze->Edges[Cur->EdgeTo[i]];
			l->Prev = NULL;
			Maze->Connected[Friend] = 1;
			Enqueue(Queue, l);
		}
	}
	while (Queue->size > 0) {
		Maze->Selected[((edge *)Queue->Head->Item)->Loc] = 1;
		Cur = GetNode(Maze, Queue->Head->V);
		for (i=0;i<4;++i) {
			Friend = Cur->Friends[i];
			if (Friend > 0) { 
				if (Maze->Connected[Friend] != 1){
					mst_link_t * l = malloc(sizeof(mst_link_t));
					l->V = Friend;
					l->Item = &Maze->Edges[Cur->EdgeTo[i]];
					l->Prev = NULL;
					Maze->Connected[Friend] = 1;
					Enqueue(Queue, l);
				}
			}
		}
		Dequeue(Queue);
	}	
	free(Queue);
	/*PyPrintEdge(Maze);*/
	return;
}
void PrimMST(circular_maze * Maze) {
	int i, Friend, Iterations;
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
			HeapEnqueue(Queue, Maze->Edges[Cur->EdgeTo[i]].Weight, &Maze->Edges[Cur->EdgeTo[i]]);
		}
	}
	
	Iterations = 0;
	while ((Item = (edge *) HeadDequeue(Queue)) != NULL) {
		Maze->Selected[Item->Loc] = 1;
		Cur = GetNode(Maze, Item->Out);
		for (i=0;i<4;++i) {
			Friend = Cur->Friends[i];
			if (Friend > 0) { 
				if (Maze->Connected[Friend] != 1){
					Maze->Edges[Cur->EdgeTo[i]].Out = Friend;
					Maze->Connected[Friend] = 1;
					HeapEnqueue(Queue, Maze->Edges[Cur->EdgeTo[i]].Weight, &Maze->Edges[Cur->EdgeTo[i]]);
				}
			}
		}
		++Iterations;
	}	
	End = omp_get_wtime();
	printf("%f:%d\n", End-Start, Iterations);
	free(Queue->Links);
	free(Queue);
	return;
}
void PrimMST1(circular_maze * Maze) {
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
void PrimMST2(circular_maze * Maze) {
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
	while ((Item = (edge *) HeadDequeue1(Queue)) != NULL) {
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
void PrimMST3(circular_maze * Maze) {
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
			if (Friend > 0) { 
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
    PrimMST3(&Maze);
    /*cleanup*/
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

    return 0;
}