
/*
Nome: Keegan Carlo
Cognome: Falcao
Matricola: 0001042769
Classe: A
Indirizzo email: keegancarlo.falcao@studio.unibo.it
*/

/*IMPORTANT
APPROACH EXPLAINED:

    As the input is a matrix, I structured the nodes as a matrix of coords,
    instead of a simple array of nodes.
    I made this choice to easily access to bidimensional array (matrix) 
    i.e nodes, previous, distance, visited values such as
        → barray[coord.x][coord.y]

    Further we need to print the path as coordinate format and such approach grant to be pretty useful

*/
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <assert.h>
#include <limits.h>
#include <math.h>

#define INPUT 4
#define MAX_SIZE 250
#define MIN_SIZE 5
#define MAX_ADJAX 4

#define CHECK(I) printf(" ~Checkpoint %d\n", I);
#define EVAL(V) printf("\nValues is → %d ", V);
#define NEWLINE() printf("\n");

/*STRUCTS*/
typedef struct
{
    int **barray;
    int n;
    int m;
} Matrix;

typedef struct
{
    int x;
    int y;
} Coord;

typedef struct
{
    Coord source;
    Coord destination;
    long int weight;
} Edge;

typedef struct
{
    Coord coord;
    Edge *edges;
} Node;

typedef struct
{
    Node **nodes;
    int n;
    int m;
} Graph;

/*METHODS DECLARATIONS*/
/*read*/
int read_file(FILE **file, const char *input_file_path);
int read_int(FILE *file, int *value);
int read_long(FILE *file, long int *value);
/*Reads an input matrix structure
[INPUT]: 
-Input file
-output matrix (pointer to pointer)
-number of columns
-number of rows
[RETURN]: control bit 
*/
int read_matrix(FILE *file, int **mat, int m, int n);

/*initializations*/
Matrix *init_input_matrix(FILE *file, long int *cf, long int *ch);
/*Creates a graph from a matrix of nodes
[INPUT]: 
-Matrix to extrapolate its values and its shape
-height difference needed to compute the weight from a node to a successive one
[RETURN]: Graph pointer
*/
Graph *create_graph(Matrix *H, long int ch);
void init_node(Node *node, int i, int j);
/*Initialize each node with its adjax nodes (MAX 4), the graph becomes is represented .. 
as matrix where each node has an adjiancy list represented by the edges (naturally edges ..
..overlap but, it is not big deal)
[INPUT]: 
-i, j → current node coordinates
-graph → to check the adjiax bounds 
-barray → values to compute the weight
-ch → height difference to compute weight
[RETURN]: \
*/
void init_adjax(int i, int j, Graph *g, int **barray, long int ch);
int insert_edge(Edge *edge, int *c, long int weight, Coord destination, Coord source);
/*Initializes distance as a bidimensional array to Infinite ~To later get the minimum~
[INPUT]: 
-bidimensional array which contains best weight from source
-Graph to get distance array shape
[RETURN]: Graph pointer
*/
void init_dist(long int** barr, Graph* g);

/*creation*/
Coord **create_matrix_coord(Graph* g);
long int **create_matrix_longint(Graph* g);
int **create_matrix_int(Graph* g);

/*destruction*/
void destroy_matrix(Matrix *mat);
void destroy_graph(Graph *g);
void destroy_node(Node *node);
void destroy_matrix_int(int **mat, int rows);
void destroy_matrix_longint(long int **mat, int rows);
void destroy_matrix_coord(Coord **mat, int rows);

/*weight*/

/*Computes the weight from the given formula
-current value
-succesive value
-height difference multiplier
*/
long int compute_weight(int curr, int succ, long int ch);

/*coord utils*/
Coord create_coord(int x, int y);
void display_coord(Coord c);
/*Compares two coord and returns 0 if are not equal 1 otherwise*/
int is_equal_coord(Coord c1, Coord c2);

/*chekers*/
int check_bounds(int x, int y, int boundx, int boundy);
int check_param(int argc, const char *prog);

/*dijkstra helpers*/
/*Calcultes the minum coordinate excluding the visited ones
[INPUT]:
-dist → distance bidimensional array to calculate the minimum
-visited → visited coordinates
-graph to get the shape 
[RETURN]: minimum coordinate
*/
Coord min_dist_coord_sofar(long int **dist, int **visited, Graph *g);
/*Recursive algorithm to print the coordinates of the shortest path
It recursively gets the previous coordinate starting from the destination coord
*/
void shortest_path(Coord src, Coord dest, Graph *g, long int **dist, Coord **prev);
void print_result(Coord src, Coord dest, long int **dist, Coord **prev, Graph *g, long int cf);

/*dijkstra procedure*/
/*Explanation shown in the procedure*/
int dijkstra(Graph *g, Coord source, Coord destination, Coord **prev_barray, long int **distance_barray, long int weight_supplier);

/*METHODS / PROCEDURES*/
/*handle input file*/
int read_int(FILE *file, int *value) {
    return fscanf(file, "%d ", value);
}

int read_long(FILE *file, long int *value) {
    return fscanf(file, "%li ", value);
}

int read_matrix(FILE *file, int **mat, int m, int n) {
    int i, j;
    int nline = 0;

    for (i = 0; i < n; i++){
        for (j = 0; j < m; j++){
            nline += read_int(file, &mat[i][j]);
        }   
    }
    assert(nline == n*m);
    return EXIT_SUCCESS;
}

int check_param(int argc, const char *prog) {
    /*Ensures test input and program name are the only arguments passed*/
    if (argc != 2) {
        fprintf(stderr, "Invocare il programma con: %s test?.in\n", prog);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int read_file(FILE **file, const char *input_file_path) {
    if (strcmp(input_file_path, "-")) {
        *file = fopen(input_file_path, "r");
        if (*file == NULL) {
            fprintf(stderr, "Can not open %s\n", input_file_path);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

/*free*/
void destroy_matrix_int(int **mat, int rows) {
    int i;
    for (i = 0; i < rows; i++)
        free(mat[i]);
    free(mat);
}

void destroy_matrix(Matrix *mat) {
    destroy_matrix_int(mat->barray, mat->n);
    free(mat);
}

void destroy_node(Node *node) {
    free(node);
}

void destroy_graph(Graph *g) {
    int i, j;
    for (i = 0; i < g->n; i++){
        for (j = 0; j < g->m; j++){
            free(g->nodes[i][j].edges);     
        }
        free(g->nodes[i]);
    }
    free(g->nodes);
    free(g);
}

void destroy_matrix_longint(long int **mat, int rows) {
    int i;
    for (i = 0; i < rows; i++)
        free(mat[i]);
    free(mat);
}

void destroy_matrix_coord(Coord **mat, int rows) {
    int i;
    for (i = 0; i < rows; i++)
        free(mat[i]);
    free(mat);
}

/*coord*/
Coord create_coord(int x, int y) {
    Coord *coord;
    Coord ret;
    coord = (Coord*)malloc(sizeof(Coord));
    coord->x = x;
    coord->y = y;
    ret = *coord;
    free(coord);
    return ret;
}

void display_coord(Coord c) {
    printf("%d %d ", c.x, c.y);
}

int is_equal_coord(Coord c1, Coord c2) {
    return c1.x == c2.x && c2.y == c1.y;
}

long int compute_weight(int curr, int succ, long int ch) {
    return ch * (curr - succ) * (curr - succ);
}

/*INITIALIZATION*/
Matrix *init_input_matrix(FILE *file, long int *cf, long int *ch) {
    int i;
    int nline;
    int m, n;
    Matrix *H;
    int **barray;

    assert(file != NULL);
    nline = 0;
    nline += read_long(file, cf);
    nline += read_long(file, ch);
    nline += read_int(file, &n);
    nline += read_int(file, &m);

    assert(nline == INPUT);
    assert(n <= MAX_SIZE && n >= MIN_SIZE);
    assert(m <= MAX_SIZE && m >= MIN_SIZE);
    
    /*allocate memory to the matrix, each row in presented as a single array pointer, so a matrix is a double pointer*/
    barray = (int**)malloc((n) * sizeof(int*));
    assert(barray != NULL);
    
    for (i = 0; i < n; i++) {
        barray[i] = (int*)malloc(m * sizeof(int));
        assert(barray[i] != NULL);  
    }
    /*Allocate memory for matrix*/
    H = malloc(sizeof(Matrix));
    assert(H != NULL);

    /*assignament*/
    read_matrix(file, barray, m, n);
    H->barray = barray;
    H->n = n;
    H->m = m;

    return H;
}

void init_node(Node *node, int i, int j) {
    node->coord = create_coord(i, j);
    node->edges = calloc(MAX_ADJAX, sizeof(Edge));
    assert(node->edges != NULL);    
}

int insert_edge(Edge *edge, int *c, long int weight, Coord destination, Coord source) {
    edge[*c].source = source;
    edge[*c].destination = destination;
    edge[*c].weight = weight;
    *c += 1;        
    return EXIT_SUCCESS;
}

void init_adjax(int i, int j, Graph *g, int **barray, long int ch) {
    int startx, starty;
    int adjiax;
    int c;
    Coord dest;

    c = 0;
    
    for(startx = i - 1; startx <= i + 1; startx++) {
        for (starty = j - 1; starty <= j + 1; starty++) {
            if (check_bounds(startx, starty, g->n, g->m) && (startx == i || starty == j) && !(startx == i && starty == j)) {
                adjiax = barray[startx][starty];
                dest = create_coord(startx, starty);
                insert_edge(g->nodes[i][j].edges, &c, compute_weight(barray[i][j], adjiax, ch), dest, g->nodes[i][j].coord);  
            }
        }
    }
}

void init_dist(long int** barr, Graph* g) {
    int i, j;
    for (i = 0; i < g->n; i++)
        for (j = 0; j < g->m; j++)
            barr[i][j] = LONG_MAX;
}

/*outbound check*/
int check_bounds(int x, int y, int boundx, int boundy) {
    return x >= 0 && x < boundx && y >= 0 && y < boundy;
}

/*memalloc*/
Graph *create_graph(Matrix *H, long int ch) {
    int **barray;
    int i, j, v;
    Graph *g;
    
    barray = H->barray;
    g = malloc(sizeof(Graph));
    g->m = H->m;
    g->n = H->n;
    g->nodes = malloc(g->n * sizeof(Node*));
    assert(g->nodes != NULL);
    for (v = 0; v < g->n; v++) {
        g->nodes[v] = (Node*)malloc(g->m * sizeof(Node));
        assert(g->nodes[v] != NULL);
    }
    /*Allocate memory for matrix*/

    for (i = 0; i < H->n; i++){
        for (j = 0; j < H->m; j++){
            init_node(&g->nodes[i][j], i, j);
            init_adjax(i, j, g, barray, ch);  
        }
    }
    return g;
}

Coord **create_matrix_coord(Graph* g) {
    Coord **mat;
    int v, j;

    mat = (Coord**)malloc(g->n * sizeof(Coord*));
    assert(mat != NULL);
    for (v = 0; v < g->n; v++) {
        mat[v] = (Coord*)malloc(g->m * sizeof(Coord));
        assert(mat[v] != NULL);
        for(j = 0; j < g->m; j++) {
            mat[v][j] = create_coord(-1, -1);
        }
    }
    return mat;
}

long int **create_matrix_longint(Graph* g) {
    long int **mat;
    int v, i;

    mat = (long int**)malloc(g->n * sizeof(long int*));
    assert(mat != NULL);
    for (v = 0; v < g->n; v++) {
        mat[v] = (long int*)malloc(g->m * sizeof(long int));
        assert(mat[v] != NULL);
        for (i = 0; i < g->m; i++){
            mat[v][i] = 0;
        }
    }
    return mat;
}

int **create_matrix_int(Graph* g) {
    int **mat;
    int v, i;

    mat = (int**)malloc(g->n * sizeof(int*));
    assert(mat != NULL);
    for (v = 0; v < g->n; v++) {
        mat[v] = (int*)malloc(g->m * sizeof(int));
        assert(mat[v] != NULL);
        for (i = 0; i < g->m; i++){
            mat[v][i] = 0;
        }   
    }
    return mat;
}

/*dijkstra*/
Coord min_dist_coord_sofar(long int **dist, int **visited, Graph *g) {
    int i, j;
    long int min_val;
    Coord min_edge;

    min_val = LONG_MAX;
    for (i = 0; i < g->n; i++) {
        for (j = 0; j < g->m; j++) {
            int exist = visited[i][j];
            if(!exist) {
                if(dist[i][j] < min_val) {
                    min_val = dist[i][j];
                    min_edge = create_coord(i, j);
                }
            }
        }   
    }
    return min_edge; 
}

void shortest_path(Coord src, Coord dest, Graph *g, long int **dist, Coord **prev) {
    if (is_equal_coord(dest, src))
        return;
    dest = prev[g->nodes[dest.x][dest.y].coord.x][g->nodes[dest.x][dest.y].coord.y];
    shortest_path(src, dest, g, dist, prev);
    display_coord(dest);
    NEWLINE()
    return;
}

void print_result(Coord src, Coord dest, long int **dist, Coord **prev, Graph *g, long int cf) {
    long int totalsum;

    shortest_path(src, dest, g, dist, prev);
    display_coord(dest);
    NEWLINE()
    totalsum = dist[dest.x][dest.y] + cf;
    /*mark ending*/
    display_coord(create_coord(-1, -1));
    NEWLINE()
    printf("%li",totalsum);
}

int dijkstra(Graph *g, Coord src, Coord dest, Coord **prev, long int **dist, long int weight_supp) {
    Coord curr;
    int **visited;
    Edge *edges;
    int c;
    long int accumulate_dist;
    Coord adjax_node;

    /*Initial coordinate is 0, 0*/
    curr = src;
    visited = create_matrix_int(g);
    init_dist(dist, g);
    dist[curr.x][curr.y] = 0;
    prev[curr.x][curr.y] = curr;

    /*Time complexity O(V) where V is the number of vertexes because is the worst case it may..
    ..lear all nodes weight until it found the destination*/
    /*Overall Time Complexity O(V) * [O(V) ~minimum algorithm~] → O(V^2)*/
    while(!is_equal_coord(curr, dest)) {
        edges = g->nodes[curr.x][curr.y].edges;
        visited[curr.x][curr.y] = 1;
    
        for (c = 0; c < MAX_ADJAX; c++)
        {
            adjax_node = edges[c].destination;
            /*checks for the adjax only if not visited */
            if (!visited[adjax_node.x][adjax_node.y]) {
                /*accumlated distance to the adjax_{c} node  */
                accumulate_dist = edges[c].weight + (weight_supp) + dist[curr.x][curr.y];
                /*checks whether the current unvisited node weight is less than..
                ..the previous so it chooses the best one for that particular node
                */
                if (accumulate_dist < dist[adjax_node.x][adjax_node.y]) {
                    prev[adjax_node.x][adjax_node.y] = curr;
                    dist[adjax_node.x][adjax_node.y] = accumulate_dist;
                }
            }    
        }
        /*Minimum is found from a traditional array so..*/
        /*Time complexity is O(E) where E are the number of vertexes, in this case n*m */
        curr = min_dist_coord_sofar(dist, visited, g);
        if (is_equal_coord(curr, create_coord(-1, -1))) {
            fprintf(stderr, "Coordinata minima non adeguata (-1, -1)");
            return EXIT_FAILURE;                
        }     
    }
    destroy_matrix_int(visited, g->n);
    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) 
{ 
    FILE *filein = stdin;
    /*input matrix, taken from file*/
    Matrix *H;
    long int diff_cost;
    long int fixed_cost;
    /*main graph structered as a adjency list, each node has its coord, and all its edges*/
    Graph *g;
    Coord src, dest;
    /*previous nodes (represented by coords) matrix*/
    Coord **prev;
    /*distance cost of all nodes*/
    long int **dist;
    /*assures that functions do not stumble upon errors*/
    int result;

    result = check_param(argc, argv[0]);
    if (result) return EXIT_FAILURE;

    result = read_file(&filein, argv[1]);
    if (result) return EXIT_FAILURE;
    
    H = init_input_matrix(filein, &fixed_cost, &diff_cost);
    g = create_graph(H, diff_cost);
    
    src = create_coord(0, 0);
    dest = create_coord(g->n-1, g->m-1);
    dist = create_matrix_longint(g);
    prev = create_matrix_coord(g);
    
    result = dijkstra(g, src, dest, prev, dist, fixed_cost);
    if (result) return EXIT_FAILURE;

    print_result(src, dest, dist, prev, g, fixed_cost);
    NEWLINE()

    destroy_matrix_longint(dist, g->n);
    destroy_matrix_coord(prev, g->n);
    destroy_graph(g);
    destroy_matrix(H);
    
    if (filein != stdin) fclose(filein);
	return EXIT_SUCCESS; 
}