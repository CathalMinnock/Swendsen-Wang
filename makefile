OBJECTS = main.o grid.o 
CC = scorep mpicc
CFLAGS = -Wall -lm -L '/home/support/apps/cports/rhel-6.x86_64/gnu/papi/5.6.0/lib' -lmpi
EXEC = swprog

make: $(OBJECTS)
	$(CC) -o $(EXEC) $(OBJECTS) $(CFLAGS)

grid.o: grid.c grid.h
	$(CC) -c grid.c $(CFLAGS)
	
parallel_main.o: main.c grid.h
	$(CC) -c main.c $(CFLAGS)
	
clean:
	rm -f $(OBJECTS) $(EXEC)
