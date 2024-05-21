/* C code to implement shortest path from src to destination given height differences */

/*Author: Keegan Carlo Falcao
Institution: ALMA Mater Studiorium Bologna (Unibo)
Difficulty: medium
*/

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <assert.h>
#include <limits.h>

#define INPUT 4
#define MAX_SIZE 250
#define MIN_SIZE 5
#define MAX_ADJAX 4

#define CHECK(I) printf("Checkpoint %d\n", I);
#define EVAL(V) printf("\nValues is → %d: ", V);
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
int read_matrix(FILE *file, int **mat, int m, int n);

/*initializations*/
Matrix *init_input_matrix(FILE *file, long int *cf, long int *ch);
Graph *create_graph(Matrix *H, long int ch);
int init_node(Node *node, int i, int j);
int init_adjax(int i, int j, Graph *g, int **barray, long int ch);
int insert_edge(Edge *edge, int *c, long int weight, Coord destination, Coord source);
int init_dist(long int** barr, Graph* g);

/*creation*/
Coord **create_matrix_coord(Graph* g);
long int **create_matrix_longint(Graph* g);
int **create_matrix_int(Graph* g);

/*destruction*/
int destroy_matrix(Matrix *mat);
int destroy_graph(Graph *g);
int destroy_node(Node *node);
int destroy_matrix_int(int **mat, int rows);
int destroy_matrix_longint(long int **mat, int rows);
int destroy_matrix_coord(Coord **mat, int rows);

/*weight*/
int init_weight(Edge *edges);
long int compute_weight(int curr, int succ, long int ch);

/*coord utils*/
Coord create_coord(int x, int y);
int display_coord(Coord c);
int cmp_coords(Coord c1, Coord c2);

/*chekers*/
int check_bounds(int x, int y, int boundx, int boundy);
int check_param(int argc, const char *prog);

/*dijkstra helpers*/
Coord min_dist_coord_sofar(long int **dist, int **visited, Graph *g);
int shorterst_path(Coord src, Coord dest, Graph *g, long int **dist, Coord **prev);
int print_result(Coord src, Coord dest, long int **dist, Coord **prev, Graph *g, long int cf);

/*dijkstra procedure*/
int dijkstra(Graph *g, Coord src, Coord dest, Coord **prev, long int **dist, long int weight_supp);

/*METHODS*/
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
    return 0;
}

int check_param(int argc, const char *prog) {
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
int destroy_matrix_int(int **mat, int rows) {
    int i;
    for (i = 0; i < rows; i++)
        free(mat[i]);
    free(mat);
    return 0;
}

int destroy_matrix(Matrix *mat) {
    destroy_matrix_int(mat->barray, mat->n);
    free(mat);
    return 0;
}

int destroy_node(Node *node) {
    free(node);
    return 0;
}

int destroy_graph(Graph *g) {
    int i, j;
    for (i = 0; i < g->n; i++){
        for (j = 0; j < g->m; j++){
            free(g->nodes[i][j].edges);     
        }
        free(g->nodes[i]);
    }
    free(g->nodes);
    free(g);
    
    return 0;
}

int destroy_matrix_longint(long int **mat, int rows) {
    int i;
    for (i = 0; i < rows; i++)
        free(mat[i]);
    free(mat);
    return 0;
}

int destroy_matrix_coord(Coord **mat, int rows) {
    int i;
    for (i = 0; i < rows; i++)
        free(mat[i]);
    free(mat);
    return 0;
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

int display_coord(Coord c) {
    printf("%d %d ", c.x, c.y);
    return 0;
}

int cmp_coords(Coord c1, Coord c2) {
    return c1.x == c2.x && c2.y == c1.y;
}

/*weight*/
int init_weight(Edge *edges) {
    int i;
    for (i = 0; i < MAX_ADJAX; i++){
        edges[i].weight = INT_MAX;  
    }
    return 0;
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

int init_node(Node *node, int i, int j) {
    node->coord = create_coord(i, j);
    node->edges = calloc(MAX_ADJAX, sizeof(Edge));
    assert(node->edges != NULL);
    
    /*initialize edges weight*/
    init_weight(node->edges);
    return 0;
}

int insert_edge(Edge *edge, int *c, long int weight, Coord destination, Coord source) {
    
    edge[*c].source = source;
    edge[*c].destination = destination;
    edge[*c].weight = weight;
    *c += 1;        
    return 0;
}

int init_adjax(int i, int j, Graph *g, int **barray, long int ch) {
    int startx, starty;
    int adjiax;
    int c;
    Coord dest;

    c = 0;
    if (cmp_coords(create_coord(i, j), create_coord(g->n-1, g->m-1)) == 1) {
        return 0;
    }

    for(startx = i - 1; startx <= i + 1; startx++) {
        for (starty = j - 1; starty <= j + 1; starty++) {
            if (check_bounds(startx, starty, g->n, g->m)) {
                if ((startx == i || starty == j) && !(startx == i && starty == j)) {
                    adjiax = barray[startx][starty];
                    dest = create_coord(startx, starty);
                    insert_edge(g->nodes[i][j].edges, &c, compute_weight(barray[i][j], adjiax, ch), dest, g->nodes[i][j].coord);  
                }
            }
        }
    }
    return 0;
}

int init_dist(long int** barr, Graph* g) {
    int i, j;
    for (i = 0; i < g->n; i++)
        for (j = 0; j < g->m; j++)
            barr[i][j] = INT_MAX;
    return 0;
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
    long int min_val = LONG_MAX;
    Coord min_edge;
    for (i = 0; i < g->n; i++)
    {
        for (j = 0; j < g->m; j++)
        {
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

int shorterst_path(Coord src, Coord dest, Graph *g, long int **dist, Coord **prev) {
    if (cmp_coords(dest, src)) {
        return 0;
    }
    dest = prev[g->nodes[dest.x][dest.y].coord.x][g->nodes[dest.x][dest.y].coord.y];
    shorterst_path(src, dest, g, dist, prev);
    display_coord(dest);
    NEWLINE()
    return 0;
}

int print_result(Coord src, Coord dest, long int **dist, Coord **prev, Graph *g, long int cf) {
    long int totalsum;

    shorterst_path(src, dest, g, dist, prev);
    display_coord(dest);
    NEWLINE()
    totalsum = dist[dest.x][dest.y] + cf;
    /*mark ending*/
    display_coord(create_coord(-1, -1));
    NEWLINE()
    printf("%li",totalsum);
    return 0;
}

int dijkstra(Graph *g, Coord src, Coord dest, Coord **prev, long int **dist, long int weight_supp) {

    Coord curr;
    int **visited;
    Edge *edges;
    int c;
    long int accumulate_dist;
    Coord adjax_node;

    curr = src;
    visited = create_matrix_int(g);
    
    init_dist(dist, g);
    dist[curr.x][curr.y] = 0;
    prev[curr.x][curr.y] = curr;

    /*clears all nodes weight until it found the destination*/
    while(!cmp_coords(curr, dest)) {
        edges = g->nodes[curr.x][curr.y].edges;
        visited[curr.x][curr.y] = 1;
    
        for (c = 0; c < MAX_ADJAX; c++)
        {
            adjax_node = edges[c].destination;
            /*check for the adjiax onli if not visited*/
            if (!visited[adjax_node.x][adjax_node.y]) {
                accumulate_dist = edges[c].weight + (weight_supp) + dist[curr.x][curr.y];
                /*checks whether the current unvisited node weight is less than 
                the previous so it chooses the best one for that particular node
                */
                if (accumulate_dist < dist[adjax_node.x][adjax_node.y]) {
                    prev[adjax_node.x][adjax_node.y] = curr;
                    dist[adjax_node.x][adjax_node.y] = accumulate_dist;
                }
            }    
        }
        
        curr = min_dist_coord_sofar(dist, visited, g);
        if (cmp_coords(curr, create_coord(-1, -1))) {
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

    result = print_result(src, dest, dist, prev, g, fixed_cost);
    if (result) return EXIT_FAILURE;
    NEWLINE()

    
    destroy_matrix_longint(dist, g->n);
    destroy_matrix_coord(prev, g->n);
    destroy_graph(g);
    destroy_matrix(H);
    
    if (filein != stdin) fclose(filein);
	return EXIT_SUCCESS; 
}
