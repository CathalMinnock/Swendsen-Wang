#include "grid.h"

void magnetization() {
	int i, j, k;
	double x_sum = 0, y_sum = 0;
	for(i = 0; i < local_x_size; ++i)
	for(j = 0; j < local_y_size; ++j)
	for(k = 0; k < local_z_size; ++k) {
		x_sum += x_values[ lattice[i][j][k].spin ];
		y_sum += y_values[ lattice[i][j][k].spin ];
	}
	x_sum /= x_size * y_size * z_size;
	y_sum /= x_size * y_size * z_size;
	double total_x_sum, total_y_sum;
	MPI_Reduce(&x_sum, &total_x_sum, 1, MPI_DOUBLE, MPI_SUM, root, CART_COMM);
	MPI_Reduce(&y_sum, &total_y_sum, 1, MPI_DOUBLE, MPI_SUM, root, CART_COMM);
	if(rank == root) {
		mag = sqrt(total_x_sum * total_x_sum + total_y_sum * total_y_sum);
	}
}

void sw_iterate() {
	reset_lattice();
	send_lattice_to_boundaries();
	setup_bond_configuration();
	reduce_local_labels();
	
	total_finished = false;
	while(total_finished == false) {
		local_finished = true;
		
		send_lattice_to_boundaries_x();
		reduce_boundary_labels_x();
		send_boundaries_to_lattice_x();
		
		send_lattice_to_boundaries_y();
		reduce_boundary_labels_y();
		send_boundaries_to_lattice_y();
		
		send_lattice_to_boundaries_z();
		reduce_boundary_labels_z();
		send_boundaries_to_lattice_z();
		
		reduce_local_labels();
		MPI_Allreduce(&local_finished, &total_finished, 1, MPI_C_BOOL, MPI_LAND, CART_COMM);
	}
	
	generate_rand_cluster_spins();
	flip_clusters();
}

void reset_lattice() {
	int i, j, k;
	for(i = 0; i < local_x_size; ++i)
	for(j = 0; j < local_y_size; ++j)
	for(k = 0; k < local_z_size; ++k) {
		lattice[i][j][k].label = (i + s[0]) * (y_size * z_size) + (j + s[1]) * z_size + (k + s[2]);
		lattice[i][j][k].bond_x = false;
		lattice[i][j][k].bond_y = false;
		lattice[i][j][k].bond_z = false;
	}
}

void setup_bond_configuration() {
	int i, j, k;
	for(i = 0; i < local_x_size; ++i)
	for(j = 0; j < local_y_size; ++j)
	for(k = 0; k < local_z_size; ++k) {
		if(lattice[i][j][k].spin == lattice[i+1][j][k].spin && (double) rand() / (double) RAND_MAX < prob)
			lattice[i][j][k].bond_x = true;
		if(lattice[i][j][k].spin == lattice[i][j+1][k].spin && (double) rand() / (double) RAND_MAX < prob)
			lattice[i][j][k].bond_y = true;
		if(k == local_z_size - 1) {
			if(lattice[i][j][k].spin == front[i][j].spin && (double) rand() / (double) RAND_MAX < prob)
				lattice[i][j][k].bond_z = true;
		}
		else if(lattice[i][j][k].spin == lattice[i][j][k+1].spin && (double) rand() / (double) RAND_MAX < prob)
			lattice[i][j][k].bond_z = true;
	}
}

// Needs to know bonds, updates labels
// reduce cluster labels using just bonds between local points, NOT across boundaries
void reduce_local_labels() {
	int i, j, k;
	bool finished = false;
	while(!finished) {
		finished = true;
		// X direction bonds
		for(i = 0; i < local_x_size - 1; ++i)
		for(j = 0; j < local_y_size; ++j)
		for(k = 0; k < local_z_size; ++k) {
			if(lattice[i][j][k].bond_x) {
				if(lattice[i][j][k].label > lattice[i+1][j][k].label) {
					lattice[i][j][k].label = lattice[i+1][j][k].label;
					finished = false;
				}
				else if(lattice[i+1][j][k].label > lattice[i][j][k].label) {
					lattice[i+1][j][k].label = lattice[i][j][k].label;
					finished = false;
				}
			}
		}
		// Y direction bonds
		for(i = 0; i < local_x_size; ++i)
		for(j = 0; j < local_y_size - 1; ++j)
		for(k = 0; k < local_z_size; ++k) {
			if(lattice[i][j][k].bond_y) {
				if(lattice[i][j][k].label > lattice[i][j+1][k].label) {
					lattice[i][j][k].label = lattice[i][j+1][k].label;
					finished = false;
				}
				else if(lattice[i][j+1][k].label > lattice[i][j][k].label) {
					lattice[i][j+1][k].label = lattice[i][j][k].label;
					finished = false;
				}
			}
		}
		// Z direction bonds
		for(i = 0; i < local_x_size; ++i)
		for(j = 0; j < local_y_size; ++j)
		for(k = 0; k < local_z_size - 1; ++k) {
			if(lattice[i][j][k].bond_z) {
				if(lattice[i][j][k].label > lattice[i][j][k+1].label) {
					lattice[i][j][k].label = lattice[i][j][k+1].label;
					finished = false;
				}
				else if(lattice[i][j][k+1].label > lattice[i][j][k].label) {
					lattice[i][j][k+1].label = lattice[i][j][k].label;
					finished = false;
				}
			}
		}
	}
}
void send_lattice_to_boundaries() {
	MPI_Request req[6];
	MPI_Isend(&lattice[0][0][0], 1, TB_LATTICE_PLANE, top_p, 0, CART_COMM, &req[0]); // Send to Top
	MPI_Irecv(&bottom[0][0], 1, TB_BOUNDARY_PLANE, bottom_p, 0, CART_COMM, &req[1]); // Recv from Bottom
	MPI_Isend(&lattice[0][0][0], 1, LR_LATTICE_PLANE, left_p, 1, CART_COMM, &req[2]); // Send to Left
	MPI_Irecv(&right[0][0], 1, LR_BOUNDARY_PLANE, right_p, 1, CART_COMM, &req[3]); // Recv from Right
	MPI_Isend(&lattice[0][0][0], 1, FB_LATTICE_PLANE, back_p, 2, CART_COMM, &req[4]); // Send to Back
	MPI_Irecv(&front[0][0], 1, FB_BOUNDARY_PLANE, front_p, 2, CART_COMM, &req[5]); // Recv from Front
	MPI_Waitall(6, req, MPI_STATUSES_IGNORE);
}
void send_lattice_to_boundaries_x() {
	MPI_Request req[2];
	MPI_Isend(&lattice[0][0][0], 1, TB_LATTICE_PLANE, top_p, 0, CART_COMM, &req[0]); // Send to Top
	MPI_Irecv(&bottom[0][0], 1, TB_BOUNDARY_PLANE, bottom_p, 0, CART_COMM, &req[1]); // Recv from Bottom
	MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
}
void send_lattice_to_boundaries_y() {
	MPI_Request req[2];
	MPI_Isend(&lattice[0][0][0], 1, LR_LATTICE_PLANE, left_p, 1, CART_COMM, &req[0]); // Send to Left
	MPI_Irecv(&right[0][0], 1, LR_BOUNDARY_PLANE, right_p, 1, CART_COMM, &req[1]); // Recv from Right
	MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
}
void send_lattice_to_boundaries_z() {
	MPI_Request req[2];
	MPI_Isend(&lattice[0][0][0], 1, FB_LATTICE_PLANE, back_p, 2, CART_COMM, &req[0]); // Send to Back
	MPI_Irecv(&front[0][0], 1, FB_BOUNDARY_PLANE, front_p, 2, CART_COMM, &req[1]); // Recv from Front
	MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
}

void send_boundaries_to_lattice() {
	MPI_Request req[6];
	MPI_Isend(&bottom[0][0], 1, TB_BOUNDARY_PLANE, bottom_p, 0, CART_COMM, &req[0]); // Send to Bottom
	MPI_Irecv(&lattice[0][0][0], 1, TB_LATTICE_PLANE, top_p, 0, CART_COMM, &req[1]); // Recv from Top
	MPI_Isend(&right[0][0], 1, LR_BOUNDARY_PLANE, right_p, 1, CART_COMM, &req[2]); // Send to Right
	MPI_Irecv(&lattice[0][0][0], 1, LR_LATTICE_PLANE, left_p, 1, CART_COMM, &req[3]); // Recv from Left
	MPI_Isend(&front[0][0], 1, FB_BOUNDARY_PLANE, front_p, 2, CART_COMM, &req[4]); // Send to Front
	MPI_Irecv(&lattice[0][0][0], 1, FB_LATTICE_PLANE, back_p, 2, CART_COMM, &req[5]); // Recv from Back
	MPI_Waitall(6, req, MPI_STATUSES_IGNORE);
}
void send_boundaries_to_lattice_x() {
	MPI_Request req[2];
	MPI_Isend(&bottom[0][0], 1, TB_BOUNDARY_PLANE, bottom_p, 0, CART_COMM, &req[0]); // Send to Bottom
	MPI_Irecv(&lattice[0][0][0], 1, TB_LATTICE_PLANE, top_p, 0, CART_COMM, &req[1]); // Recv from Top
	MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
}
void send_boundaries_to_lattice_y() {
	MPI_Request req[2];
	MPI_Isend(&right[0][0], 1, LR_BOUNDARY_PLANE, right_p, 1, CART_COMM, &req[0]); // Send to Right
	MPI_Irecv(&lattice[0][0][0], 1, LR_LATTICE_PLANE, left_p, 1, CART_COMM, &req[1]); // Recv from Left
	MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
}
void send_boundaries_to_lattice_z() {
	MPI_Request req[2];
	MPI_Isend(&front[0][0], 1, FB_BOUNDARY_PLANE, front_p, 2, CART_COMM, &req[0]); // Send to Front
	MPI_Irecv(&lattice[0][0][0], 1, FB_LATTICE_PLANE, back_p, 2, CART_COMM, &req[1]); // Recv from Back
	MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
}

void reduce_boundary_labels_x() {
	int j, k;
	// Bottom boundary
	for(j = 0; j < local_y_size; ++j)
	for(k = 0; k < local_z_size; ++k) {
		if(lattice[local_x_size - 1][j][k].bond_x) {
			if(bottom[j][k].label > lattice[local_x_size - 1][j][k].label) {
				bottom[j][k].label = lattice[local_x_size - 1][j][k].label;
				local_finished = false;
			}
			else if(lattice[local_x_size - 1][j][k].label > bottom[j][k].label) {
				lattice[local_x_size - 1][j][k].label = bottom[j][k].label;
				local_finished = false;
			}
		}
	}
}
void reduce_boundary_labels_y() {
	int i, k;
	// Right boundary
	for(i = 0; i < local_x_size; ++i)
	for(k = 0; k < local_z_size; ++k) {
		if(lattice[i][local_y_size - 1][k].bond_y) {
			if(right[i][k].label > lattice[i][local_y_size - 1][k].label) {
				right[i][k].label = lattice[i][local_y_size - 1][k].label;
				local_finished = false;
			}
			else if(lattice[i][local_y_size - 1][k].label > right[i][k].label) {
				lattice[i][local_y_size - 1][k].label = right[i][k].label;
				local_finished = false;
			}
		}
	}
}
void reduce_boundary_labels_z() {
	int i, j;
	// Front boundary
	for(i = 0; i < local_x_size; ++i)
	for(j = 0; j < local_y_size; ++j) {
		if(lattice[i][j][local_z_size - 1].bond_z) {
			if(front[i][j].label > lattice[i][j][local_z_size - 1].label) {
				front[i][j].label = lattice[i][j][local_z_size - 1].label;
				local_finished = false;
			}
			else if(lattice[i][j][local_z_size - 1].label > front[i][j].label) {
				lattice[i][j][local_z_size - 1].label = front[i][j].label;
				local_finished = false;
			}
		}
	}
}

void generate_rand_cluster_spins() {
	int i;
	for(i = 0; i < recvcounts[rank]; ++i)
		local_rand_spins[i] = rand() % q;
	MPI_Allgatherv(local_rand_spins, recvcounts[rank], MPI_UNSIGNED_CHAR, rand_spins, recvcounts, displacements, MPI_UNSIGNED_CHAR, CART_COMM);
}
void flip_clusters() {
	int i, j, k;
	for(i = 0; i < local_x_size; ++i)
	for(j = 0; j < local_y_size; ++j)
	for(k = 0; k < local_z_size; ++k) {
		lattice[i][j][k].spin = rand_spins[ lattice[i][j][k].label ];
	}
}
// Creates the 3D lattice, and the 3 boundary planes (for now I believe the other 3 are not needed)
void create_local_grid()
{	
	bottom = create_2D_array(local_y_size, local_z_size);		
	right = create_2D_array(local_x_size, local_z_size);	
	front = create_2D_array(local_x_size, local_y_size);	
	lattice = create_grid();
}

point **create_2D_array(int rows, int cols) {
	point *array_data = malloc(rows * cols * sizeof(point) );
	point **array = malloc(rows * sizeof(point *) );
	int i;
	for(i = 0; i < rows; ++i)
		array[i] = &(array_data[i * cols]);
	return array;
}
void free_2D_array(point **array, int rows) {
	free(array[0]);
	free(array);
}

// Builds in left / right and top / bottom boundary. Front and back boundary will need to be otherwise handled.
// Initialises to random values.
point ***create_grid() {
	lattice = malloc( (local_x_size + 1) * sizeof(point **) );
	point *lattice_data = malloc(local_x_size * local_y_size * local_z_size * sizeof(point) );
	int i, j, k;
	for(i = 0; i < local_x_size; ++i) {
		lattice[i] = malloc( (local_y_size + 1) * sizeof(point *) );
		for(j = 0; j < local_y_size; ++j) { 
			lattice[i][j] = &lattice_data[i * (local_y_size * local_z_size) + j * local_z_size];
			for(k = 0; k < local_z_size; ++k) { 
				lattice[i][j][k].spin = rand() % q; 
			} 
		}
		lattice[i][local_y_size] = right[i];
	}
	lattice[local_x_size] = bottom;
	reset_lattice(); // init the labels and bonds
	return lattice;
}
void free_grid() {
	free(lattice[0][0]);
	int i;
	for(i = 0; i < local_x_size; ++i)
		free(lattice[i]); 
	free(lattice);
}

// Figures out cart_dims[] for the processor
// Figures out cart_coordinates[] for the processor
// Figures out the start and end position for each processor in each direction: s[0], s[1], s[2], e[0], e[1], e[2]
void decomp3d(int *s, int *e)
{
	int a, b, c;
	a = cbrt(nprocs);
	while(a > 0) {
		if(nprocs % a == 0) {
			b = nprocs / a;
			c = sqrt(b);
			while(c > 0) {
				if(b % c == 0) {
					cart_dims[0] = a;
					cart_dims[1] = c;
					cart_dims[2] = b / c;
					break;
				}
				--c;
			}
			break;
		}
		--a;
	}
	cart_coordinates[0] = rank / ( cart_dims[1] * cart_dims[2]);
	cart_coordinates[1] = (rank / cart_dims[2]) % cart_dims[1];
	cart_coordinates[2] = rank % cart_dims[2];
	
	int x[3]; // Some of the processors will receive x values in each dimension, others will receive x+1
	x[0] = x_size / cart_dims[0];
	x[1] = y_size / cart_dims[1];
	x[2] = z_size / cart_dims[2];
	int remainder[3]; // The number to receive x+1 in each dimension is given by this remainder
	remainder[0] = x_size % cart_dims[0]; 
	remainder[1] = y_size % cart_dims[1];
	remainder[2] = z_size % cart_dims[2];
	int i;
	for(i = 0; i < 3; ++i) {
		// The ones to receive x+1 values are set up to be at the start
		if(cart_coordinates[i] < remainder[i]) {
			s[i] = (x[i] + 1) * cart_coordinates[i];
			e[i] = s[i] + x[i] + 1;
		}
		// Everything after them receives x values
		else {
			s[i] = (x[i] + 1) * remainder[i] + x[i] * (cart_coordinates[i] - remainder[i]);
			e[i] = s[i] + x[i];
		}
	}
	local_x_size = e[0] - s[0]; 
	local_y_size = e[1] - s[1]; 
	local_z_size = e[2] - s[2]; 
}

void print_grid() {
	int i,j,k;
	for(k=0; k < local_z_size; ++k) {
		printf("========== Rank %i: k = %i ==========\n", rank, k);
		for(i=0; i < local_x_size; ++i){
			for(j=0; j < local_y_size; ++j)
				printf("%i ", lattice[i][j][k].spin);
			printf("\n");
		}
		printf("\n");
	}
}
void print_grid_labels() {
	int i,j,k;
	for(k=0; k < local_z_size; ++k) {
		printf("========== Rank %i: k = %i ==========\n", rank, k);
		for(i=0; i < local_x_size; ++i){
			for(j=0; j < local_y_size; ++j)
				printf("%i ", lattice[i][j][k].label);
			printf("\n");
		}
		printf("\n");
	}
}
void print_grid_bond_x() {
	int i,j,k;
	for(k=0; k < local_z_size; ++k) {
		printf("========== Rank %i: k = %i ==========\n", rank, k);
		for(i=0; i < local_x_size; ++i){
			for(j=0; j < local_y_size; ++j)
				printf("%i ", lattice[i][j][k].bond_x);
			printf("\n");
		}
		printf("\n");
	}
}
void print_2D_array(point **array, int rows, int cols) {
	int i, j;
	for(i=0; i < rows; ++i){
		for(j=0; j < cols; ++j)
			printf("%i ", array[i][j].label);
		printf("\n");
	}
}
