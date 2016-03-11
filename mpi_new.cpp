#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common_mpi.h"
#include <mpi.h>
#include "libmpi.h"
#include <iostream>

//
//  benchmarking program
//

#define density 0.0005
#define mass	0.01
#define cutoff	0.01
#define min_r	(cutoff/100)
#define dt	0.0005

using namespace std;

// direction Directions[9] = {{-1,-1},{-1,0},{-1,1},{0,-1},{0,0},{0,1},{1,-1},{1,0},{1,1}};


int main( int argc, char **argv )
{   
    int navg,nabsavg=0; 
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg, rdmin;
    int rnavg;
	
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" ); 
        printf( "-no turns off all correctness checks and particle output\n");   
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    int n_proc, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;      

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous(6, MPI_DOUBLE, &PARTICLE);
    MPI_Type_commit(&PARTICLE);

    int particle_per_proc = (n + n_proc - 1) / n_proc;
    int *partition_offsets = (int*) malloc((n_proc+1) * sizeof(int));
    for(int i = 0; i < n_proc+1; i++)
        partition_offsets[i] = min(i*particle_per_proc, n);
    
    int *partition_sizes = (int*) malloc(n_proc * sizeof(int));
    for(int i = 0; i < n_proc; i++)
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
   
    int nlocal = partition_sizes[rank];
    particle_t *local = (particle_t*) malloc(nlocal * sizeof(particle_t));

    set_size( n );
  //  if(rank == 0)
  
      init_particles( n, particles );
    
      double gridSize = sqrt(density * n);
      int binNumPerEdge = ceil(gridSize / cutoff);
      int binNumber = binNumPerEdge * binNumPerEdge;
      bin_t* bins = new bin_t[binNumber];
      int* binIdxForParticle = (int*) malloc(n *sizeof(int));
    
    MPI_Scatterv(particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD);
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        davg = 0.0;
	dmin = 1.0;
        
        MPI_Allgatherv(local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD);        

        for (int i = 0; i < binNumber; i++)
        {
            bins[i].particlesInBin.clear();
        }
      //  cout << "moew1" << endl;
        for(int i = 0; i < n; i++)
        {
            int bin_x = particles[i].x / cutoff;
            int bin_y = particles[i].y / cutoff;
            int binIndex = bin_x + bin_y * binNumPerEdge;
            bins[binIndex].particlesInBin.push_back(particles+i);
            binIdxForParticle[i] = binIndex;
        }
        
	//
        //  compute all forces
        //
      //  cout << "meow2" << endl;
        for( int i = 0; i < nlocal; i++)
        {
            apply_force_per_particle(local[i], bins, binIdxForParticle[i+partition_offsets[rank]], binNumPerEdge, &dmin, &davg, &navg);
        }
        
       // cout << "meow3" << endl;
        //
        //  move particles
        //

        for( int i = 0; i < nlocal; i++ ) 
            move( local[i] );
  
        if( find_option( argc, argv, "-no" ) == -1 ) 
        {
          if(fsave && (step%SAVEFREQ) == 0)
            save(fsave, n, particles);

          MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
	  //
          //  compute statistical data
          //
          if(rank == 0){
            if (rnavg) { 
              absavg += rdavg/rnavg;
              nabsavg++;
            }
	    if (rdmin < absmin) absmin = rdmin; 
	  }	
        
        }
    }
    simulation_time = read_timer( ) - simulation_time;
    if(rank == 0) { 
      printf( "n = %d,n_proc = %d, simulation time = %g seconds", n,n_proc, simulation_time);

      if( find_option( argc, argv, "-no" ) == -1 )
      {
        if (nabsavg) absavg /= nabsavg;
    // 
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
      printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
      if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
      if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
      }
      printf("\n");
   
    //
    // Printing summary data
    //
      if( fsum)
          fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );

    free( particles );
    delete[] bins;
    free(partition_offsets);
    free(partition_sizes);
    free(binIdxForParticle);
    

    if( fsave )
        fclose( fsave );
    MPI_Finalize();
    return 0;
}
