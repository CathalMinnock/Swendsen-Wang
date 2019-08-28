#include "grid.h"

void init(int argc, char *argv[])
{
	srand(time(NULL));
	
	int option_index = 0;
	x_size = 8; y_size = 8; z_size = 8; 
	q = 2; beta = 0.4; filename = "apple.txt";
	samples = 1; steps_between_samples = 1;
	while((option_index = getopt(argc, argv, "x:y:z:q:b:s:f:" )) != -1) {
		switch(option_index) {
			case 'x':
				x_size = atoi(optarg);
				break;
			case 'y':
				y_size = atoi(optarg);
				break;
			case 'z':
				z_size = atoi(optarg);
				break;
			case 'q':
				q = atoi(optarg);
				break;
			case 'b':
				beta = atof(optarg);
				break;
			case 's':
				samples = atoi(optarg);
				break;
                        case 'a':
                                steps_between_samples = atoi(optarg);
                                break;
			case 'f':
				filename = optarg;
				break;
			default:
				printf("Option incorrect\n");
				exit(1);
		}
	}
	fp = fopen(filename,"w");
	
	// Set up the Cartesian MPI Structure
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	printf("Rank %d: %d %d %d %d %lf %s\n", rank, x_size, y_size, z_size, q, beta, filename);
	decomp3d(s, e);
	int periodic_dims[3] = {1, 1, 1};
	MPI_Cart_create(MPI_COMM_WORLD, 3, cart_dims, periodic_dims, 0, &CART_COMM);
	MPI_Cart_shift(CART_COMM, 0, 1, &top_p, &bottom_p);
	MPI_Cart_shift(CART_COMM, 1, 1, &left_p, &right_p);
	MPI_Cart_shift(CART_COMM, 2, 1, &back_p, &front_p);
	root = nprocs - 1;
	
	// Make an MPI Datatype for the point struct
	MPI_Datatype point_type;
	MPI_Datatype oldtypes[3] = {MPI_INT, MPI_UNSIGNED_CHAR, MPI_C_BOOL};
	int blockcounts[3] = {1, 1, 3};
	MPI_Aint exINT, exCHAR;
	MPI_Type_extent(MPI_INT, &exINT);
	MPI_Type_extent(MPI_UNSIGNED_CHAR, &exCHAR);
	
	MPI_Aint offsets[3] = {0, exINT, exINT + exCHAR};
	MPI_Type_create_struct(3, blockcounts, offsets, oldtypes, &point_type);
	MPI_Type_commit(&point_type);
	

	// Make MPI Datatypes for the boundary & lattice planes (top/bottom, front/back, left/right)
	MPI_Type_vector(local_y_size * local_z_size, 1, 1, point_type, &TB_BOUNDARY_PLANE);
	MPI_Type_commit(&TB_BOUNDARY_PLANE);
	MPI_Type_vector(local_x_size * local_y_size, 1, 1, point_type, &FB_BOUNDARY_PLANE);
	MPI_Type_commit(&FB_BOUNDARY_PLANE);
	MPI_Type_vector(local_x_size * local_z_size, 1, 1, point_type, &LR_BOUNDARY_PLANE);
	MPI_Type_commit(&LR_BOUNDARY_PLANE);
	
	MPI_Type_vector(local_y_size * local_z_size, 1, 1, point_type, &TB_LATTICE_PLANE);
	MPI_Type_commit(&TB_LATTICE_PLANE);
	MPI_Type_vector(local_x_size, local_z_size, local_y_size * local_z_size, point_type, &LR_LATTICE_PLANE);
	MPI_Type_commit(&LR_LATTICE_PLANE);
	MPI_Type_vector(local_x_size * local_y_size, 1, local_z_size, point_type, &FB_LATTICE_PLANE);
	MPI_Type_commit(&FB_LATTICE_PLANE);
	
	
	// Create the local lattice, its boundary planes, and initialize them
	create_local_grid();
	send_lattice_to_boundaries();
	
	// Label "lookup table"
	recvcounts = malloc(nprocs * sizeof(int));
	displacements = malloc(nprocs * sizeof(int));
	int remainder = (x_size * y_size * z_size) % nprocs;
	int x = (x_size * y_size * z_size) / nprocs;
	int i;
	for(i=0; i < nprocs; ++i) {
		if(i < remainder) {
			recvcounts[i] = x+1;
			displacements[i] = i * (x + 1);
		}
		else {
			recvcounts[i] = x;
			displacements[i] = remainder * (x + 1) + (i - remainder) * x;
		}
	}
	rand_spins = malloc(x_size * y_size * z_size * sizeof(unsigned char));
	local_rand_spins = malloc(recvcounts[rank] * sizeof(unsigned char));
	
	// Pre-calculate all the sines and cosines I need beforehand
	x_values = malloc(q * sizeof(double));
	y_values = malloc(q * sizeof(double));
	for(i=0; i < q; ++i) {
		x_values[i] = sin( (2 * M_PI * i ) / (double) q);
		y_values[i] = cos( (2 * M_PI * i ) / (double) q);
	}
	prob = 1 - exp(-2 * beta);
}

void finalize() {
	fclose(fp);
	free_grid();
	free_2D_array(bottom, local_y_size);
	free_2D_array(right, local_x_size);
	free_2D_array(front, local_x_size);
	free(rand_spins);
	free(local_rand_spins);
	free(recvcounts);
	free(displacements);
	free(x_values);
	free(y_values);
	MPI_Finalize();
}

int main(int argc, char *argv[])
{
	init(argc, argv);
	int i,j;
	for(i=0; i < steps_between_samples * 10; ++i) 
		sw_iterate();
	for(i=0; i < samples; ++i) {
		for(j=0; j < steps_between_samples; ++j)
			sw_iterate();
		magnetization();
		if(rank == root)
			fprintf(fp, "%lf\n", mag);
	}
	finalize();
	return 0;
}
