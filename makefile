OBJECTS = main.o grid.o 
CC = mpicc
CFLAGS = -Wall -lm 
EXEC = swprog

make: $(OBJECTS)
	$(CC) -o $(EXEC) $(OBJECTS) $(CFLAGS)

grid.o: grid.c grid.h
	$(CC) -c grid.c $(CFLAGS)
	
parallel_main.o: main.c grid.h
	$(CC) -c main.c $(CFLAGS)
	
clean:
	rm -f $(OBJECTS) $(EXEC)
