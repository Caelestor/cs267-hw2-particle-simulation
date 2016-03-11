#
# Edison - NERSC 
#
# Intel Compilers are loaded by default; for other compilers please check the module list
#
CC = CC
MPCC = CC
OPENMP = -openmp #Note: this is the flag for Intel compilers. Change this to -fopenmp for GNU compilers. See http://www.nersc.gov/users/computational-systems/edison/programming/using-openmp/
CFLAGS = -O3
LIBS =


TARGETS = autograder openmp_new mpi_new serial_new

all:	$(TARGETS)

#serial: serial.o common.o
#	$(CC) -o $@ $(LIBS) serial.o common.o
serial_new: serial_new.o common.o world.o grid.o
	$(CC) -o $@ $(LIBS) serial_new.o common.o world.o grid.o
autograder: autograder.o common.o
	$(CC) -o $@ $(LIBS) autograder.o common.o
#pthreads: pthreads.o common.o
#	$(CC) -o $@ $(LIBS) -lpthread pthreads.o common.o
#openmp: openmp.o common.o
#	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o
#mpi: mpi.o common.o
#	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o
#serial-try: serial-try.o world.o grid.o common.o
#	$(CC) -o $@ $(LIBS) serial-try.o world.o grid.o common.o
#openmp-try: openmp-try.o world_omp.o grid.o common.o
#	$(CC) -o $@ $(LIBS) $(OPENMP) openmp-try.o world_omp.o grid.o common.o
openmp_new: openmp_new.o common.o libomp.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp_new.o common.o libomp.o
mpi_new: mpi_new.o common.o libmpi.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi_new.o common.o libmpi.o


autograder.o: autograder.cpp common.h
	$(CC) -c $(CFLAGS) autograder.cpp
#openmp.o: openmp.cpp common.h
#	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp
#serial.o: serial.cpp common.h
#	$(CC) -c $(CFLAGS) serial.cpp
serial_new.o: serial_new.cpp world.h common.h grid.h world.cpp common.cpp grid.cpp
	$(CC) -c $(CFLAGS) serial_new.cpp
#openmp-try.o: openmp-try.cpp world_omp.h common.h grid.h world_omp.cpp common.cpp grid.cpp
#	$(CC) -c $(OPENMP) $(CFLAGS) openmp-try.cpp
openmp_new.o: openmp_new.cpp libomp.h libomp.cpp common.h common.cpp
	$(CC) -c $(OPENMP) $(CFLAGS) openmp_new.cpp
mpi_new.o: mpi_new.cpp libmpi.h libmpi.cpp
	$(MPCC) -c $(CFLAGS) mpi_new.cpp

#pthreads.o: pthreads.cpp common.h
#	$(CC) -c $(CFLAGS) pthreads.cpp
#mpi.o: mpi.cpp common.h
#	$(MPCC) -c $(CFLAGS) mpi.cpp

libomp.o: libomp.cpp libomp.h
	$(CC) -c $(CFLAGS) libomp.cpp
libmpi.o: libmpi.cpp libmpi.h
	$(CC) -c $(CFLAGS) libmpi.cpp

world.o: world.cpp world.h common.h common.cpp grid.h grid.cpp
	$(CC) -c $(CFLAGS) world.cpp

#world_omp.o: world_omp.cpp world_omp.h common.h common.cpp grid.h grid.cpp
#	$(CC) -c $(CFLAGS) $(OPENMP) world_omp.cpp
	
grid.o: grid.cpp grid.h common.h common.cpp
	$(CC) -c $(CFLAGS) grid.cpp

common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp
	

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
