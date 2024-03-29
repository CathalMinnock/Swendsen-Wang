#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h> 
#include <time.h>
#include <unistd.h>

/*
 * The point struct represents all information about one site of the lattice
 * labels are used for identifying which cluster the site is in
 * spin is a number in {0,1,...,q-1}
 * bond_x/y/z represent whether or not the site is bonded with its nearest neighbour in each of these (positive) directions
 */
typedef struct point {
	int label; 
	unsigned char spin; 
	bool bond_x; 
	bool bond_y;
	bool bond_z;
} point;

point ***create_grid(); // allocates memory for a lattice, initialises spins to random values, connects up the right & bottom boundary planes to the lattice (NOT the front plane)
point **create_2D_array(int rows, int cols); // used to make the 2D boundary planes
void create_local_grid(); // uses the above functions to make the local sub-lattice and its 3 boundary planes (front, right, bottom)
void decomp3d(int *s, int *e); // For each processor, gets cart_dims[], cart_coordinates[], s[], e[] (the start and end positions of the sub-lattice within the total lattice) 
void free_grid();
void free_2D_array(point **array, int rows);

void sw_iterate(); // performs one iteration of the Swendsen-Wang algorithm, using this set of functions
void reset_lattice(); // resets the labels of each site to their initial values, and sets all the bonds to false
void setup_bond_configuration(); // sets up the bond variables for every site
void reduce_local_labels(); // makes every cluster in the local lattice have the same label

// makes the sites in the x-direction boundary plane its neighbouring sites in the local lattice have matching labels
void reduce_boundary_labels_x();
void reduce_boundary_labels_y(); 
void reduce_boundary_labels_z();

// sends the data on the positive x, y, and z boundaries of the local lattice to the corresponding boundaries planesof the neighbouring
void send_lattice_to_boundaries();
void send_lattice_to_boundaries_x();
void send_lattice_to_boundaries_y();
void send_lattice_to_boundaries_z();

// sends the data in the x, y, and z boundary planes to the corresponding points in the neighbouring local lattice
void send_boundaries_to_lattice();
void send_boundaries_to_lattice_x();
void send_boundaries_to_lattice_y();
void send_boundaries_to_lattice_z();

void generate_rand_cluster_spins(); // for each possible cluster label, generate a random spin for sites with that cluster label to flip their spin to.
void flip_clusters(); // flip each point of the local lattice to the spin in the rand_spins array corresponding to its label
void magnetization();

// Print functions are just for debugging and testing
void print_grid();
void print_grid_labels();
void print_2D_array(point **array, int rows, int cols);
void print_grid_bond_x();
void print_grid_bond_y();
void print_grid_bond_z();

point ***lattice; // The local 3D Lattice
point **bottom, **right, **front; // Boundary planes
int x_size, y_size, z_size; // size in each dimension of the total lattice
int local_x_size, local_y_size, local_z_size; // size in each dimension of the processor's local sub-lattice

int nprocs; // number of processors in total
int cart_dims[3]; // number of processers in each direction 
int cart_coordinates[3]; // location of processor in cartesian communicator
int root; // the root rank
int rank; // rank in MPI_COMM_WORLD
int top_p, bottom_p, left_p, right_p, back_p, front_p; // ranks of neighbouring processors
int s[3], e[3]; // start and end location of the processor's data in the full lattice
MPI_Comm CART_COMM; // cartesian communicator
MPI_Datatype point_type; // MPI Datatype for the point struct
MPI_Datatype TB_LATTICE_PLANE, FB_LATTICE_PLANE, LR_LATTICE_PLANE, TB_BOUNDARY_PLANE, FB_BOUNDARY_PLANE, LR_BOUNDARY_PLANE; // MPI Datatypes for lattice planes and boundary planes

unsigned char q; // number of spin states
double beta; // inverse temperature
double prob; // probability of a bond forming between two sites
int samples; // the number of samples to take
int steps_between_samples; // the number of Swendsen-Wang iterations to perform between each sample 
int thermalize_steps; // the number of Swendsen-Wang iterations to perform before starting to take samples
double *x_values, *y_values; // arrays of pre-calculated sines and cosines of all the relevant angles

unsigned char *rand_spins; // array of N spins: each cluster chooses one of them to flip to based on its cluster label
unsigned char *local_rand_spins; // each processor's section of the full array of N spins, to be broadcasted to everyone else using Allgatherv()
int *recvcounts; // array of how many spins each processor sends
int *displacements; // array of locations in rand_spins where each processor's spins start
bool local_finished; // if the local lattice has not changed in the latest labeling iteration
bool global_finished; // if local_finished is true for every processor after an iteration
 
double mag; // the root processor stores the result of magnetization() here
char *filename; // the file to print magnetizations to
FILE *fp;
