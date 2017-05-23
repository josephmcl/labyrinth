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

typedef struct Link mst_link_t; //TODO: Header?
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
void Enqueue (mst_queue_t * Queue, mst_link_t * Item) { //TODO: make safe
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
void Dequeue (mst_queue_t * Queue) { //TODO: make safe
	if (Queue->size > 0) {	
		mst_link_t * Temp = Queue->Head;
		if (Queue->Head->Prev != NULL) {
			Queue->Head = Queue->Head->Prev;
			free(Temp);
			Queue->size -= 1;

		}
		else if (Queue->Tail != NULL) { //NOT SAFE AT ALL!
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
    if (h->Length + 1 >= h->size) {
        h->size = h->size ? h->size * 2 : 4;
        h->Links = realloc(h->Links, h->size * sizeof(mst_link_t));
    }
    int i = h->Length + 1;
    int j = i / 2;
    while (i > 1 && h->Links[j].V > Priority) {
        h->Links[i] = h->Links[j];
        i = j;
        j = j / 2;
    }
    h->Links[i].V = Priority;
    h->Links[i].Item = Item;
    h->Length++;
}
void * HeadDequeue (mst_heap_t *h) {
    int i, j, k;
    if (!h->Length) {
        return NULL;
    }
    void * Item = h->Links[1].Item;
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

void PrintEdge(edge Edge) {
	printf("%f`%f`%f`%f", Edge.A->X, Edge.A->Y, Edge.B->X, Edge.B->Y);
}
void PrintEdges(edge * Edges, int Limit) {
	int i;
	for (i=0;i<Limit;++i)
		printf("((%f,%f),(%f,%f))", Edges[i].A->X, Edges[i].A->Y, Edges[i].B->X, Edges[i].B->Y);
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
	int i, j;
	for(i=1;i<=Maze->Radius;i++) {
		float Angle = (180.0*acos((2.0*(i*i)-1.0)/(2.0*(i*i))))/pi;
		if (Angle*2.0f < Maze->Delta) {
			Maze->CurrentRingSize *= 2;
			Maze->Delta = Angle;
		}
		Maze->Angles[i-1] = Maze->Delta;
		Maze->Points[i-1] = malloc(Maze->CurrentRingSize * sizeof(int)); 
		int Top = Maze->_Size + Maze->CurrentRingSize;
		//#pragma omp parallel for
		for (j=0;j<Maze->CurrentRingSize;j++)
			Maze->Points[i-1][j] = Maze->_Size + j;
		Maze->_Size = Top;
		Maze->Sections[i-1] = Maze->CurrentRingSize;
	}
	Maze->Connected = malloc(Maze->_Size * sizeof(int));
	//#pragma omp parallel for
	for (i=0;i<Maze->_Size;i++)
		Maze->Connected[i] = 0; /* Initialize node represenation to unconnected */
	Maze->Nodes = malloc(Maze->_Size * sizeof(node *));
	int l, k = 0;
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
	k = 0;
	long NumEdges = 0;
	NumEdges = 0;
	for (i=0;i<Maze->Radius-1;++i)
		NumEdges += Maze->Sections[i] * 2;
	NumEdges += Maze->Sections[i];
	Maze->Edges = malloc(NumEdges * sizeof(edge));
	for (j=0;j<Maze->Sections[0];++j) {
		Maze->Nodes[0][j].Friends[3] = -1;
	}
	#pragma omp parallel for
	for (i=0;i<(Maze->Radius-1);++i) {
		k = PriorSum(Maze->Sections, i);
		//printf("S: %d, %d, %d\n", Maze->Sections[i], i, k);
		for (j=0;j<Maze->Sections[i];++j) {
			/* set "theta" edges. (The ones that make rings.) */
			Maze->Edges[k + j].Weight = (k + j) % THETA_WEIGHT;
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
			Maze->Edges[k + j].Weight = (k + j) % THETA_WEIGHT;
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
		
		Maze->Edges[k + j].Weight = (k + j) % THETA_WEIGHT;
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
	Maze->Selected = malloc(sizeof(int) * NumEdges);
	//#pragma omp parallel for
	for (i=0;i<NumEdges;i++) 
		Maze->Selected[i] = 0;
	return;
}

int * PrimMST(circular_maze * Maze) {

	int i, r, b, Friend;
	int a[4] = {2,1,3,0};
	node * Cur;
	mst_queue_t * Queue = malloc(sizeof(mst_queue_t)); 
	Queue->Empty = 1; 
	Queue->size = 0;
	Cur = GetNode(Maze, Maze->_Size/2);
	b = rand() % 4;
	for (r=0;r<4;++r)
		a[r] = (b+r)%4;
	for (r=0;r<4;++r) {
		i = a[r];
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
		//PrintLinks(Queue->Head);
		Maze->Selected[((edge *)Queue->Head->Item)->Loc] = 1;
		Cur = GetNode(Maze, Queue->Head->V);
		//printf("Node %d w/ friends {", Queue->Head->V);
		b = rand() % 4;
			for (r=0;r<4;++r)
			a[r] = (b+r)%4;
		for (r=0;r<4;++r) {
			i = a[r];
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

	for(i=Maze->NumEdges-1;i>0;--i){
		if(Maze->Selected[i]) {
			PrintEdge(Maze->Edges[i]);
			printf(",");
		}
	}
	PrintEdge(Maze->Edges[0]);
}
int * PrimMST1(circular_maze * Maze) {
	int i, r, b, Friend;
	node * Cur;
	edge * Item;
	mst_heap_t * Queue = malloc(sizeof(mst_heap_t)); 
	Queue->Links = malloc(sizeof(mst_link_t));
	Queue->Length = 0;
	Queue->size = 0;
	Cur = GetNode(Maze, 0);
	Maze->Connected[0] = 1;
	Maze->Selected[0] = 1;
	b = rand() % 4;
	for (i=0;i<4;++i) {
		Friend = Cur->Friends[i];
		if (Friend != -1 && Maze->Connected[Friend] != 1){
			Maze->Edges[Cur->EdgeTo[i]].Out = Friend;
			HeapEnqueue(Queue, Maze->Edges[Cur->EdgeTo[i]].Weight, &Maze->Edges[Cur->EdgeTo[i]]);
		}
	}
	while ((Item = (edge *) HeadDequeue(Queue)) != NULL) {
		//PrintLinks(Queue->Head);
		//printf("%d\n", Item->Loc);
		Maze->Selected[Item->Loc] = 1;
		Cur = GetNode(Maze, Item->Out);
		//printf("Node %d w/ friends {", Queue->Head->V);

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
	}	
	free(Queue->Links);
	free(Queue);
	for(i=Maze->NumEdges-1;i>0;--i){
		if(Maze->Selected[i]) {
			PrintEdge(Maze->Edges[i]);
			printf(",");
		}
	}
	PrintEdge(Maze->Edges[0]);

}
int * PrimMSTP(circular_maze * Maze) {
	
}

int main(int argc, char** argv) {
    if (argc != 2) {
        printf("Usage: %s <radius>\n",argv[0]);
        return 1;
    }
    // TODO: >1 err checking 
    circular_maze Maze;
    Maze.CurrentRingSize 
    				= 6; /* Should the user be able to define this? */
    Maze.Radius 	= atoi(argv[1])-1;
    Maze.Delta 		= 360.0f/6.0f;
    Maze.Points 	= malloc(Maze.Radius * sizeof(int *)); 
    Maze.Sections 	= malloc(Maze.Radius * sizeof(int *)); 
    Maze.Angles 	= malloc(Maze.Radius * sizeof(float *)); 
    Maze._Size 		= 0;
    InitalizeNodes(&Maze);
    InitalizeEdges(&Maze);
    //PrintEdges(Maze.Edges, Maze.NumEdges);
    //PrimMST(&Maze);
    //PrimMST1(&Maze);
    /*cleanup*/
    int i;
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

    return 0;
}