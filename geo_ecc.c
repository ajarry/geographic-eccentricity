/*
 * This program computes the geographic eccentricity and various other measures
 * of 2-dimensional networks under various communication models. See
 * http://arxiv.org/abs/1403.3007
 *
 * Usage: [repeat x]   repeat experiment x times
 *        [L x]        simulates 4x² nodes scattered in a 4x * x rectangle
 *        [error x]    adds a position error with Gaussian distribution and
 *           standard deviation x
 *        [truncate x] removes links that are apparently longer than x
 *        [rand x]     uses the random graph model with link probability x
 *        [sinr x y]   uses the SINR model with min range x and max range y
 *        [exp x]      uses the expontential link model with average range x
 *        [file x]     appends compound results to the end of file x
 *
 * This program invokes "voronoi", which is Steven Fortune's C program
 * for computing 2-dimensional Voronoi diagrams, available at his
 * webpage http://ect.bell-labs.com/who/sjf/
 * The voronoi program is invoked with the command:
 *  system("./voronoi -t <points.tmp >delaunay.tmp");
 * The two files points.tmp and delaunay.tmp are thus created.
 *
 * This program uses rand()/RAND_MAX to simulate a uniform random number
 * distribution in [0,1]. It also uses the Box-Muller-Marsaglia transform to
 * generate a Gaussian distribution from it. This may lead to the Neave anomaly 
 * depending on the C implementation (I didn't test it).
 *
 * Compile this program with '-lm' to link in the math library.
 *
 * Author: Aubin Jarry
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int size=2500;      /* Number of nodes. */

double* x; /* x coordinates of the nodes */
double* y; /* y coordinates of the nodes */

#define RAND 0
#define SINR 1
#define EXP 2
int model = SINR;   /* communication model */
double p = 0.004;   /* link probability for the RAND model */
double r_min = 1;   /* minimum range for the SINR model */
double r_max = 5; /* maximum range for the SINR model */
double r_avg = 1.7; /* average range for the EXP model before truncation. */

double r_err = 0;    /* standard error */ 
double r_tr = 10000; /* links with apparent length > r_tr are discarded. */

#define MAX_AVG_DEGREE 100 /* The maximum average degree, needed for reserving
			    * (*adj_buf) memory. */

int* matrix;   /* half distance matrix of the communication graph */
int* dist_buf; /* buffer of currently computed distances,
		* (compute_distances) */
int** adj_buf; /* adjacency lists of the communication graph,
		* (compute_distances) */
int*  degree;  /* degree of nodes:
		  - in the communication graph, (compute_distances)
		  - in the Delaunay triangulation, (init_voronoi embed_graph) */

int*  delaunay_raw; /* raw input from Fortune's program (init_voronoi) */
int** delaunay_dat; /* adjacency lists of the Delaunay triangulation,
		    * (init_voronoi embed_graph) */

int diameter;      /* communication graph diameter, >=size means infinite */
int delaunay_loc;  /* locality of the Delaunay triangulation */
int geo_ecc;       /* geographic eccentricity */
int embedding_loc; /* locality of the link embedding */

clock_t t1,t2; /* measures time taken by the compute_distance function */

/****************************************************
 * Section 1: Initialization                        *
 ****************************************************/

/*
 * Safe malloc function, used only in 'init'.
 */
void* safe_malloc(int s) {
  void* ptr = malloc(s);
  if (ptr == NULL) {
    printf("Out of memory.");
    exit(-1);
  }
  return ptr;
}
/*
 * Allocates memory and initializes the random number generator.
 */
void init() {
  srand(time(NULL));
  matrix = (int*) safe_malloc(sizeof(int)*size*(size-1)/2);
  dist_buf = (int*) safe_malloc(sizeof(int)*size*(size-1));
  adj_buf = (int**) safe_malloc(sizeof(int*)*size);
  adj_buf[0] = (int*) safe_malloc(sizeof(int)*size*MAX_AVG_DEGREE);
  x = (double*) safe_malloc(sizeof(double)*size);
  y = (double*) safe_malloc(sizeof(double)*size);
  /* A Delaunay triangulation features at most (2n-4) triangles with 3 vertices
   * each, so we need 3*(2n-4) slots in delaunay_raw. */
  delaunay_raw = (int*) safe_malloc(sizeof(int)*(6*size-12));
  degree = (int*) safe_malloc(sizeof(int)*size);
  delaunay_dat = (int**) safe_malloc(sizeof(int*)*size);
  /* There are at most (3n-6) edges in the Delaunay triangulation, so the
   * adjacency lists will take up to 2*(3n-6) slots. However, we need an
   * additional slot per vertex in our implementation (see function
   * init_voronoi), which makes (7n-12) required slots. */
  delaunay_dat[0] = (int*) safe_malloc(sizeof(int)*(7*size-12)); // 2*(3n-6)
}

/****************************************************
 * Section 2: Communication graph/matrix functions. *
 ****************************************************/

/*
 * States that all distances in the graph are infinite (inifinity = size).
 */
void init_matrix() {
  int c;
  for (c=size*(size-1)/2-1; c>=0; c--) matrix[c] = size;
}
/*
 * Index of the matrix cell that correspons to the distance (u,v). The contract
 * is that u and v have to be distinct.
 */
#define c(u,v) (u<v? v*(v-1)/2+u : u*(u-1)/2+v)
/*
 * The matrix cell that corresponds to the distance (u,v) if u < v.
 */
#define m(u,v) matrix[v*(v-1)/2 + u]
/*
 * The matrix cell that corresponds to the distance (u,v) if u != v.
 */
#define mm(u,v) matrix[c(u,v)]
/*
 * Computes and stores the distances between u and w, by taking into account
 * the already computed distance (u,v) and all the edges vw. This code
 * belongs exclusively to 'compute_distances' and it appears twice there.
 */
#define ADD_HOP(u,v) {\
  int i;\
  int* list = adj_buf[u];\
  for (i=degree[u]-1; i>=0; i--) {\
    w = list[i];\
    if (v!=w && mm(v,w) > d) {\
      mm(v,w) = d;\
      *(ptr++) = v;\
      *(ptr++) = w;\
    }\
  }\
}
/*
 * Computes the distances between nodes (hop count) in the communication graph.
 * This is the time-critical function, which is the focus of any serious
 * code optimization.
 * The current implementation is in O(n*m) as it should be. Normally, our
 * communications graphs are scarce with m ~ 10n.
 */
void compute_distances() {
  t1 = clock();
  int u,v,w;
  int *ptr = dist_buf;

  /* Phase 1. Compute the adjacency lists, and initialize the distance buffer
   * with nodes at distance 1.*/
  int *buf = adj_buf[0];
  int edge_count = 0;
  int* edge_ptr = adj_buf[0];
  for (u=0; u<size; u++) {
    degree[u] = 0;
    adj_buf[u] = edge_ptr;
    for (v=0; v<size; v++) {
      if (u==v) continue;
      if (mm(u,v) == 1) {
	if (u<v) {
	  *(ptr++) = u;
	  *(ptr++) = v;
	}
	*(edge_ptr++) = v;
	degree[u]++;
	edge_count++;
      }
    }
  }
  if ((edge_count - size*MAX_AVG_DEGREE) > 0) {
    printf("Too many edges in the graph, average degree %f.\n",
	   ((double)edge_count)/size);
    exit(-1);
  }

  /* Phase 2: compute the actual distances. */
  int* ptr_d = dist_buf;
  int* ptr_next = ptr;
  int d=2;
  while(ptr_d < ptr) {
    u = *(ptr_d++);
    v = *(ptr_d++);
    ADD_HOP(u,v)
    ADD_HOP(v,u)
    if (ptr_d == ptr_next) { d++; ptr_next = ptr; }
  }

  /* Phase 3: verify that the graph is connected and compute diameter. */
  /*   largest component diameter should be d-2 or d-3. */
  d = d-3;
  int c;
  for(c=size*(size-1)/2-1; c>=0; c--) {
    if (matrix[c] > d) d=matrix[c];
  }
  diameter = d;
  t2 = clock();
}
/*
 * Computes the average neighborhood size at /hop/ hops.
 */
double neighborhood(int hop) {
  int c;
  int total = size; // begins as the sum of 0-hop neighborhood
  for (c=size*(size-1)/2-1; c>=0; c--) {
    if (matrix[c] <= hop) total+=2;
  }
  return ((double)total) / size;
}

/*****************************************************
 * Section 3: Localization and network.              *
 *****************************************************/

/*
 * Scatter /size/ nodes randomly in a rectangle from (0,0) to (x_max,y_max).
 */
void init_positions(double x_max, double y_max) {
  int u;
  for (u=size-1; u>=0; u--) {
    x[u] = (x_max*rand()) / RAND_MAX;
    y[u] = (y_max*rand()) / RAND_MAX;
  }
}
/*
 * Connects the network according to the communication model and the nodes'
 * position. Edges are expressed as '1' in the distance matrix.
 */
void connect_network() {
  int u,v;
  double sig_min, sig_max;
  double dx,dy,ds;
  switch(model) {
  case SINR:
    sig_max = 1/(r_min * r_min);
    sig_min = 1/(r_max * r_max);
    break;
  }
  for (v=size-1; v>0; v--) {
    for (u=v-1; u>=0; u--) {
      dx = x[v] - x[u];
      dy = y[v] - y[u];
      ds = dx*dx + dy*dy;
      switch(model) {
      case SINR:
	p = ((1/ds) - sig_min)/(sig_max - sig_min);
	break;
      case EXP:
	p = exp(-sqrt(ds)/r_avg);
	break;
      }
      if (p >= 1 || p>0 && rand()<= p*RAND_MAX) m(u,v) = 1;
    }
  }
}
/*
 * Box-Muller-Marsaglia transform.
 */
double saved = NAN; /* saved generated random number for BMM. */
double box_muller_marsaglia() {
  double x1,x2,s;
  if (! isnan(saved)) {
    x1 = saved;
    saved = NAN;
    return x1;
  }
  do {
    x1 = (2.0 * rand())/RAND_MAX - 1;
    x2 = (2.0 * rand())/RAND_MAX - 1;
    s = x1*x1 + x2*x2;
  } while (s >= 1 || s == 0);
  s = sqrt(-2 * log(s)/s);
  saved = x2*s;
  return x1*s;
}
/*
 * Adds a Gaussian error to node positions, with avg value 0 and standard
 * deviation r_err. Removes links with apparent length > r_tr.
 */
void add_position_error() {
  double r, alpha;
  int u,v;
  double dx,dy,ds;
  if (r_err > 0) {
    for (u=size-1; u>=0; u--) {
      alpha = (M_PI * rand())/RAND_MAX;
      r = box_muller_marsaglia()*r_err;
      x[u] += r*cos(alpha);
      y[u] += r*sin(alpha);
    }
  }
  if (r_tr < 10000) {
    r = r_tr * r_tr;
    for (v=size-1; v>0; v--) {
      for (u=v-1; u>=0; u--) {
	dx = x[v] - x[u];
	dy = y[v] - y[u];
	ds = dx*dx + dy*dy;
	if (ds > r) m(u,v) = size;
      }
    }
  }
}

/*****************************************************
 * Section 4: Voronoi diagram and graph embedding.   *
 *****************************************************/

/*
 * Adds v as a neighbor of u.
 */
void addE(int u, int v) {
  int c;
  int* list = delaunay_dat[u];
  for (c=0; c<degree[u]; c++) {
    if (list[c] == v) return;
  }
  list[c] = v;
  degree[u]++;
}
/*
 * Adds the edges i->j and j->i.
 */
#define addEdge(i,j) addE(i,j);addE(j,i)

/*
 * Imports the Delaunay triangulation of the points x,y from Fortune's program,
 * reformats it as an adjacency list graph implementation, and computes the
 * Delaunay locality.
 */
void init_voronoi() {
  /* Phase 1: print the coordinates in a file. */
  FILE* file = fopen("points.tmp", "w");
  int u,v,w;
  for (u=0; u<size; u++) fprintf(file, "%f %f\n",x[u],y[u]);
  fclose(file);
 
  /* Phase 2: call Fortune's program. */
  system("./voronoi -t <points.tmp >delaunay.tmp");

  /* Phase 3: scan the delaunay triangulation from a file and reconstruct
   * the graph. */
  file = fopen("delaunay.tmp", "r");
  int * ptr= delaunay_raw;
  int c=0; /* 3*number of scanned triangles */
  while (fscanf(file, "%d",ptr++)==1) { c++; }
  fclose(file);
  for(u=size-1; u>=0; u--) degree[u] = 0;
  for (ptr--;ptr>=delaunay_raw; ptr--) { degree[*ptr]++; }
  int sum =0; /* required storage space for previous adjacency lists. */
  for (u=0; u<size; u++) {
    delaunay_dat[u] = delaunay_dat[0] + sum;
    /* degree[u] is the number of triangles in which u appears,
     * which is exactly the degree of u UNLESS U IS ON THE OUTER FACE.
     * The actual degree of u is thus between degree[u] and degree[u]+1.
     * We overestimate the total needed space by less than /size/. */
    sum += degree[u]+1; 
    degree[u] = 0;
  }
  for(c-=3; c>=0; c-=3) {
    u = delaunay_raw[c]; v= delaunay_raw[c+1]; w= delaunay_raw[c+2];
    addEdge(u,v); addEdge(v,w); addEdge(w,u);
  }

  /* Phase 4: compute the face diameter. */
  int d=0;
  for (u=size-1; u>=0; u--) for (c=degree[u]-1; c>=0; c--) {
      v = (delaunay_dat[u])[c];
      if (d < mm(u,v)) d=mm(u,v);
  }
  delaunay_loc = d;
}

/*
 * Embeds the communication graph links as line segments in the plane and
 * computes the voronoi cells that it traverses. The embedding locality
 * and the geographic eccentricity are measured this way.
 */
void embed_graph() {
  int a,b;
  int k=0, K=0;
  for (b = size-1; b>0; b--) {
    int c_b = b*(b-1)/2;
    for (a=b-1; a>= 0; a--) if (matrix[c_b+a] == 1) {
	/* point p along the line segment [ab] */
	double x_p = x[a], y_p = y[a];
	/* u is the site of the current voronoi cell, u_next is the next. */
	int u, u_next;
	for (u=a; u!=b; u=u_next) {
	  u_next = -1;
	  double x_pb = x[b]-x_p, y_pb = y[b]-y_p;
	  double t_min = 100.0; // p + t_min*pb
	  double x_pu = x[u]-x_p, y_pu = y[u]-y_p;
	  int * edge_list = delaunay_dat[u];
	  int i,v;
	  for (i = degree[u]-1; i>=0; i--) {
	    v = edge_list[i];
	    if (v == b) {
	      u_next = b;
	      t_min = 1.0;
	      break;
	    }
	    double x_uv = x[v]-x[u], y_uv = y[v]-y[u];
	    double uv_pb = x_uv*x_pb + y_uv*y_pb;
	    if (uv_pb <= 0) continue; /* not the good direction */
	    double uv_uv = x_uv*x_uv + y_uv*y_uv;
	    double uv_pu = x_uv*x_pu + y_uv*y_pu;
	    /* length is the ground to cover on (uv): (pu+0.5 uv)*uv */
	    double length = uv_pu + 0.5*uv_uv;
	    if (length <= 0) continue; /* already outside cell */
	    /* also, length = (t*pb)*uv */
	    double t = length / uv_pb;
	    if (t < t_min) { t_min = t; u_next = v; }
	  } 
	  x_p += t_min*x_pb;
	  y_p += t_min*y_pb;
	  int d = mm(u,u_next);
	  if (d > k) {
	    k = d;
	    if (d > K) K = d;
	  }
	  d = mm(a,u_next);
	  if (d > K) K = d;
	}
      }
  }
  embedding_loc = k; geo_ecc = K;
}

/*****************************************************
 * Section 5: Main.                                  *
 *****************************************************/

int main(int argc, char** argv) {
  /* Phase 1: read parameters. */
  int repeat=1;
  int c;
  double x_max=100, y_max=25;
  char* file_name = NULL;
  for (c=2; c<argc; c+=2) {
    char *cmd = argv[c-1], *val = argv[c];
    if (strcmp("repeat", cmd)==0) repeat = atoi(val);
    else if (strcmp("L", cmd)==0) {
      y_max = atof(val);
      x_max = 4*y_max;
      size = x_max * y_max;
    }
    else if (strcmp("error", cmd)==0) r_err = atof(val);
    else if (strcmp("truncate",cmd)==0) r_tr = atof(val);
    else if (strcmp("rand",cmd)==0) {
      model = RAND;
      p = atof(val);
    }
    else if (strcmp("sinr",cmd)==0) {
      model = SINR;
      r_min = atof(val);
      c++;
      r_max = (c<argc)? atof(argv[c]) : 0;
      if (r_max == 0) r_max = 1.4*r_min;
    }
    else if (strcmp("exp",cmd)==0) {
      model = EXP;
      r_avg = atof(val);
    }
    else if (strcmp("file",cmd)==0) file_name = val;
  }

  /* Phase 2: General init. */
  init();
  double* val = (double*) safe_malloc(sizeof(double)*9*repeat);
  double avg[9];
  double dev[9];
  const char* expl[] = {"D"," N_1"," k_T"," k_e", " N_e"," k_g"," N_g"," dk"," dN"};

  /* Phase 3: Run simulations. */
  int errcount = 0;
  for(c=repeat-1; c>=0; c--) {
    init_matrix();  
    init_positions(x_max,y_max);
    connect_network();
    add_position_error();
    compute_distances();

    if (diameter >= size) {
      printf("Not connected.\n");
      errcount++;
      if (errcount >= repeat) break;
      c++;
      continue;
    }

    init_voronoi();
    embed_graph();
 
    double* ptr = val+9*c;
    ptr[0] = diameter;
    ptr[1] = neighborhood(1);
    ptr[2] = delaunay_loc;
    ptr[3] = embedding_loc;
    ptr[4] = neighborhood(embedding_loc);
    ptr[5] = geo_ecc;
    ptr[6] = neighborhood(geo_ecc);
    ptr[7] = geo_ecc - embedding_loc;
    ptr[8] = ptr[6] - ptr[4];
    int i;
    for (i=0; i<7; i++) printf("%s %.2f", expl[i], ptr[i]);
    printf("\n");
  }

  /* Phase 4: Report on simulations. */
  if (repeat == 1) {
    printf("Time count: %lds %lds %lds.\n",
	   t1/CLOCKS_PER_SEC,t2/CLOCKS_PER_SEC,clock()/CLOCKS_PER_SEC);
    exit(0);
  }
  printf("%d disconnected graphs.\n", errcount);
  printf("%d connected graphs.\n", repeat-c-1);
  if (c >= 0) exit(0);
  for (c=8; c>=0; c--) {
    avg[c] = 0; 
    int i;
    for (i=9*(repeat-1)+c; i>=0; i-=9) avg[c] += val[i];
    avg[c] = avg[c]/repeat;
    dev[c] = 0;
    for (i=9*(repeat-1)+c; i>=0; i-=9) {
      double diff = val[i] - avg[c];
      dev[c] += diff*diff;
    }
    dev[c] = sqrt(dev[c]/repeat);
  }
  FILE * file = file_name ? fopen(file_name, "a") : stdout;

  if (file_name) {
    fprintf(file, "\nArgs:");
    for (c= 1; c<argc; c++) fprintf(file," %s", argv[c]);
    fprintf(file, "\n%d disconnected graphs.\n", errcount);
  }
  fprintf(file, "Average:");
  for (c=0; c<7; c++) fprintf(file, " %.2f", avg[c]);
  fprintf(file, "\nCompound: ");
  for (c=0; c<7; c++) fprintf(file, "%s %.2f %.2f", expl[c],avg[c],dev[c]);
  for (c=7; c<9; c++) fprintf(file, "%s %.2f",expl[c],dev[c]);
  fprintf(file, "\n");
  if (file_name) fclose(file);
  return 0;
}
