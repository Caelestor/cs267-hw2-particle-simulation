#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"
#include "libomp.h"
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
    int navg,nabsavg=0,numthreads; 
    double dmin, absmin=1.0,davg,absavg=0.0;
	
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

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;      

    int p = read_int(argc, argv, "-p", 1);
    omp_set_num_threads(p);

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );


    double gridSize = sqrt(density * n);
    int binNumPerEdge = ceil(gridSize / cutoff);
    int binNumber = binNumPerEdge * binNumPerEdge;
    bin_t* bins = new bin_t[binNumber];

    omp_lock_t* locks = (omp_lock_t*) malloc(binNumber * sizeof(omp_lock_t));
    for(int i = 0; i < binNumber; i++)
        omp_init_lock(locks+i);
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    #pragma omp parallel private(dmin) 
    {
    numthreads = omp_get_num_threads();
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        davg = 0.0;
	dmin = 1.0;
        
      //  #pragma omp for schedule(static)
	#pragma omp for
        for (int i = 0; i < binNumber; i++)
        {
            bins[i].particlesInBin.clear();
        }
      //  cout << "moew1" << endl;

     //   #pragma omp for schedule(static, 125)
      	#pragma omp for
	for(int i = 0; i < n; i++)
        {
            int bin_x = particles[i].x / cutoff;
            int bin_y = particles[i].y / cutoff;
            int binIndex = bin_x + bin_y * binNumPerEdge;
            omp_set_lock(locks + binIndex);
            bins[binIndex].particlesInBin.push_back(particles+i);
            omp_unset_lock(locks+binIndex);
        }
	//
        //  compute all forces
        //
      //  cout << "meow2" << endl;
        #pragma omp for reduction (+:navg) reduction(+:davg)
        for( int i = 0; i < binNumber; i++ )
        {
	    apply_force_per_bin(bins, i, binNumPerEdge, &dmin, &davg, &navg);  
        }
        
       // cout << "meow3" << endl;
        //
        //  move particles
        //
        #pragma omp for
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );
  
        if( find_option( argc, argv, "-no" ) == -1 ) 
        {
          //
          //  compute statistical data
          //
          #pragma omp master
          if (navg) { 
            absavg += davg/navg;
            nabsavg++;
          }

          #pragma omp critical
	  if (dmin < absmin) absmin = dmin; 
		
          //
          //  save if necessary
          //
          #pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
}
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

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
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );

    free( particles );
    delete[] bins;
    free(locks);

    if( fsave )
        fclose( fsave );
    
    return 0;
}
