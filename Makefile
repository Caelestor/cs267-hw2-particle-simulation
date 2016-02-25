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


TARGETS = serial pthreads openmp mpi autograder

all:	$(TARGETS)

# serial: serial.o common.o
#	$(CC) -o $@ $(LIBS) serial.o common.o
serial: serial-try.o common_new.o world.o grid.o
	$(CC) -o $@ $(LIBS) serial-try.o common_new.o world.o grid.o
autograder: autograder.o common.o
	$(CC) -o $@ $(LIBS) autograder.o common.o
# pthreads: pthreads.o common.o
#	$(CC) -o $@ $(LIBS) -lpthread pthreads.o common.o
openmp: openmp.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o
mpi: mpi.o common.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o

autograder.o: autograder.cpp common.h
	$(CC) -c $(CFLAGS) autograder.cpp
openmp.o: openmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp
# serial.o: serial.cpp common.h
#	$(CC) -c $(CFLAGS) serial.cpp
serial-try.o: serial-try.cpp common_new.h world.h
	$(CC) -c $(CFLAGS) serial-try.cpp

# pthreads.o: pthreads.cpp common.h
#	$(CC) -c $(CFLAGS) pthreads.cpp
mpi.o: mpi.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
# common.o: common.cpp common.h
#	$(CC) -c $(CFLAGS) common.cpp

world.o: world.cpp world.h
	$(CC) -c $(CFLAGS) world.cpp
	
grid.o: grid.cpp grid.h
	$(CC) -c $(CFLAGS) grid.cpp

common_new.o: common_new.cpp common_new.h
	$(CC) -c $(CFLAGS) common_new.cpp	

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
