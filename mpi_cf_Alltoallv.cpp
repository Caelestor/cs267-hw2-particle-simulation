#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "common.h"
#include "libomp.h"
#include <iostream>
#include <string>


using namespace std;
//
//  benchmarking program
//

#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt  0.0005

int main( int argc, char **argv )
{    
    int navg = 0, nabsavg=0;
    double dmin = 1.0,absmin=1.0,davg = 0.0,absavg=0.0;
    double rdavg,rdmin;
    int rnavg; 
    //
    //  process command line parameters
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;   //this is an append file: for all simulations!

    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );

    //set the number of bins.
    double gridSize = sqrt(density * n);
    int binNumPerEdge = ceil(gridSize / cutoff);
    int binNumber = binNumPerEdge * binNumPerEdge;
    bin_t* bins = new bin_t[binNumber];
    for (int i =0; i< binNumber; i++)
        bins[i].particlesInBin.clear();
    //array of int containing the numbers of the native bins
    int *native_bins;

    //define sides of proc grid
    int rowsOfProc = ceil(sqrt(n_proc));
    int colsOfProc = floor(sqrt(n_proc));
    int keep = 1;
    for (int i = rowsOfProc; i>0 and keep==1; i --) //if n_proc is prime or no number near sqrt(n_proc) divides n_proc, this is bad
    {
        if (n_proc%i == 0) { rowsOfProc = i; colsOfProc = n_proc/i; keep = 0;}
    }    
    //assign bins to processors through the arrays of corners
    int regRowsOfBinsPerProc = binNumPerEdge/rowsOfProc;
    int regColsOfBinsPerProc = binNumPerEdge/colsOfProc;
    if (rank ==0) cout << "The "<<n_proc << " Proc form: " << rowsOfProc<<" rows, "<< colsOfProc << " cols" ;
    if (rank ==0) cout << " with binNumPerEdge = " << binNumPerEdge<< ", regRowsOfBinsPerProc = " << regRowsOfBinsPerProc << ", regColsOfBinsPerProc = " << regColsOfBinsPerProc <<endl;

    //
    // Every processors computes its coordinates and bins
    //
    int localRowsOfBins;
    int localColsOfBins;
    int localNumOfBins;
    int topLeftBin;
    int x_proc;
    int y_proc;
    x_proc = rank / rowsOfProc; //proc are in column major!! (reverse as bins)
    y_proc = rank % rowsOfProc;
    if (y_proc == rowsOfProc - 1) // if last row of procs
    {
        if (x_proc == colsOfProc -1) // if last col of procs
        {
            localRowsOfBins = binNumPerEdge - regRowsOfBinsPerProc * (rowsOfProc-1);
            localColsOfBins = binNumPerEdge - regColsOfBinsPerProc * (colsOfProc-1);
            localNumOfBins = localRowsOfBins * localColsOfBins;
            native_bins = (int*) malloc( localNumOfBins * sizeof(int) );
        }
        else //last row but not last col
        {
            localRowsOfBins = binNumPerEdge - regRowsOfBinsPerProc * (rowsOfProc-1);
            localColsOfBins = regColsOfBinsPerProc;
            localNumOfBins = localRowsOfBins * localColsOfBins;
            native_bins = (int*) malloc( localNumOfBins * sizeof(int));
        }
    }
    else //we are not in the last row proc
    {
        if (rank / rowsOfProc == colsOfProc -1) // if we are in the last column
        {
            localRowsOfBins = regRowsOfBinsPerProc;
            localColsOfBins = binNumPerEdge - regColsOfBinsPerProc * (colsOfProc-1);
            localNumOfBins = localRowsOfBins * localColsOfBins;
            native_bins = (int*) malloc( localNumOfBins * sizeof(int));
        }
        else // we are in the middle, general case
        {
            localRowsOfBins = regRowsOfBinsPerProc;
            localColsOfBins = regColsOfBinsPerProc;
            localNumOfBins = localRowsOfBins * localColsOfBins;
            native_bins = (int*) malloc( localNumOfBins * sizeof(int));
        }
    }
    topLeftBin = y_proc * regRowsOfBinsPerProc * binNumPerEdge + x_proc * regColsOfBinsPerProc;
    for (int i = 0; i < localRowsOfBins; i++)
        for(int j = 0; j< localColsOfBins; j++)
            native_bins[i * localColsOfBins + j] = topLeftBin + i * binNumPerEdge + j;
    // //check that the implementation is correct:
    // char out[1024], *put = out;
    // for (int i = 0; i< localNumOfBins; i++)
    //     put += sprintf(put, "%d, ", native_bins[i]);
    // printf("I am proc %d and my bins are: %s\n",rank,out);
    // fflush(stdout);    

    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );
    // particles that belong to the processor and where we receive the particles
    particle_t *local = (particle_t*) malloc( 9 * n * sizeof(particle_t) );//this has to change 9 
    //this will point to the positions of the above that are Native
    int localNative[n]; //positions of the Native particles in local!!!!
    int nlocalNative = 0;
    //vector structure to assign processors. We copy the full particles to flatten later.
    vector< vector<particle_t*> > assignProc(n_proc, vector<particle_t*>(NULL));
    // particles that will be scattered and how (among the ones belonging to the processor)
    particle_t *particles_scatter = (particle_t*) malloc(9 * n * sizeof(particle_t));
    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );  //why do I need +1?
    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );  

    particle_t current_part;

    if( rank == 0 )
    { 
        init_particles( n, local );
        nlocalNative = n;
    }

    //
    //  simulate a number of time steps
    //
    double time_binning = 0.0;
    double time_binning2 = 0.0;
    double time_commu = 0.0;
    double time_computing = 0.0;
    double time_moving = 0.0;
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        //initialize stats
        // if (rank == 0) printf("\n------step%d--------\n",step);
        navg = 0;
        dmin = 1.0;
        davg = 0.0;
        // here we build what is going to be scattered and start the bining
        int x_bin ;
        int y_bin ;
        for (int i =0; i< binNumber; i++)
            bins[i].particlesInBin.clear();
        for (int p = 0; p< n_proc; p++)
            assignProc[p].clear();
        // cout << "Proc "<< rank << " has cleaned the bins and assignProc"<<endl;

        // 
        // re-count the Native parts and put the pointers to them.
        // 
        if ( (step +1)%SAVEFREQ == 0 and rank == 0)
            time_binning = read_timer();
        int nlocalNativeNew = 0;
        for (int i = 0; i < nlocalNative; i++)  //I need to also send the particles that are in the boundary!! => duplicate!
        {
            x_bin = floor(local[i].x/cutoff);
            y_bin = floor(local[i].y/cutoff);
            x_proc = x_bin/regColsOfBinsPerProc;
            if (x_proc > colsOfProc -1) x_proc = colsOfProc -1;
            y_proc = y_bin/regRowsOfBinsPerProc;
            if (y_proc > rowsOfProc -1) y_proc = rowsOfProc -1;
            // std::cout << "Proc "<<rank<<", part " << i << ": x_bin = "<< x_bin << ", y_bin = "<<y_bin << " assigned to proc: ";
            // cout.flush();

            // 
            // create the new native. Don't append the particle to its own processor!
            // 
            if ( x_proc*rowsOfProc + y_proc == rank)
            {
                localNative[nlocalNativeNew] = i; //now this is an int that tells the offset!!
                nlocalNativeNew += 1;
                // cout << "Native, ";
            }
            // else
            // {
            //     std::cout << x_proc*rowsOfProc + y_proc << ", ";
            // }
            // cout.flush();

            // 
            // push_back to all bin and proc that need particle i
            // 
            bins[x_bin + binNumPerEdge * y_bin].particlesInBin.push_back(local + i); //its bin
            assignProc[x_proc*rowsOfProc + y_proc].push_back(local + i); //its own proc
            if (x_bin%regColsOfBinsPerProc == 0 and x_bin > 0)
            {
                assignProc[(x_proc-1)*rowsOfProc + y_proc].push_back(local + i); //proc to the left
                // std::cout <<(x_proc-1)*rowsOfProc + y_proc << ", ";
                if (y_bin % regRowsOfBinsPerProc == 0 and y_bin > 0)
                {
                    assignProc[(x_proc-1) * rowsOfProc + y_proc -1].push_back(local + i); //proc to the left up
                    // std::cout << (x_proc-1) * rowsOfProc + y_proc -1 << ", ";
                }
                if (y_bin % regRowsOfBinsPerProc == regRowsOfBinsPerProc -1  and y_proc < rowsOfProc -1 )
                {
                    assignProc[(x_proc-1) * rowsOfProc + y_proc +1].push_back(local + i); //proc to the left down
                    // std::cout << (x_proc-1) * rowsOfProc + y_proc +1 << ", ";
                }
            }
            if (x_bin%regColsOfBinsPerProc == regColsOfBinsPerProc -1 and x_proc < colsOfProc -1)
            {
                assignProc[((x_proc +1)*rowsOfProc + y_proc)].push_back(local + i); // proc to the right
                // std::cout << (x_proc + 1) * rowsOfProc + y_proc  << ", ";
                if ( y_bin % regRowsOfBinsPerProc == 0 and y_bin > 0)
                {
                    assignProc[(x_proc+1) * rowsOfProc + y_proc -1].push_back(local + i); // right up
                    // std::cout << (x_proc+1) * rowsOfProc + y_proc  - 1 << ", ";
                }
                if ( y_bin % regRowsOfBinsPerProc ==  regRowsOfBinsPerProc -1 and y_proc < rowsOfProc -1)
                {
                    assignProc[(x_proc+1) * rowsOfProc + y_proc +1].push_back(local + i); // right low
                    // std::cout << (x_proc+1) * rowsOfProc + y_proc +1 << ", ";
                }
            }
            if (y_bin%regRowsOfBinsPerProc == 0 and y_bin > 0)
            {
                assignProc[x_proc * rowsOfProc + y_proc -1].push_back(local + i) ; //up
                // std::cout << x_proc * rowsOfProc + y_proc -1 << ", ";
            }
            if (y_bin % regRowsOfBinsPerProc == regRowsOfBinsPerProc -1 and y_proc < rowsOfProc -1)
            {
                assignProc[x_proc * rowsOfProc + y_proc + 1].push_back(local + i); //low
                // std::cout << x_proc * rowsOfProc + y_proc +1 << ", ";
            }
            // std::cout << std::endl;
        }
        // cout.flush();

        // 
        //eliminate the particles already in the processor so they are not sent again to himself!
        // 
        assignProc[rank].clear();
        // 
        //flatten the vector assignProc and keep partition sizes/offsets
        // 
        // cout << "The partition in Proc "<<rank<<" is: ";
        for (int proc = 0; proc < n_proc; proc++)
        {
            if (proc == 0)
                partition_offsets[proc] = 0;
            else
                partition_offsets[proc] = partition_offsets[proc -1] + partition_sizes[proc -1];
            for (int part = 0; part < assignProc[proc].size(); part ++)
                particles_scatter[partition_offsets[proc] + part ] = *(assignProc[proc][part]);
            partition_sizes[proc] = assignProc[proc].size();
            // cout << partition_sizes[proc] << ", ";
        } 
        if ( (step +1)%SAVEFREQ == 0 and rank == 0)
        {
            time_binning = read_timer() - time_binning;
            time_commu = read_timer();
        }
        // cout << "--> "<< nlocalNativeNew <<" localNatives" <<endl;
        // cout.flush();

        // //let's check that the particles are allright
        // cout << endl << "The particles gathered have for x:"<< endl;
        // for (int proc = 0; proc < n_proc; proc++)
        // {
        //     cout << "Proc " << proc <<" : ";
        //     for (int part = 0; part < assignProc[proc].size(); part ++)
        //         cout << particles_scatter[partition_offsets[proc] + part ].x << ", ";
        //     cout << endl;
        // }       

        //
        //  allocate storage for local partition
        //
        int nlocalsReceived[n_proc];
        for (int p = 0; p < n_proc; p++) nlocalsReceived[p] = 0;
        int local_offsets[n_proc];
        for (int p = 0; p < n_proc; p++) local_offsets[p] = 0;
        int ntotalReceived = 0;

        MPI_Alltoall( partition_sizes  , 1, MPI_INT, nlocalsReceived          , 1, MPI_INT, MPI_COMM_WORLD );
        // 
        //compute nlocal
        // 
        // char out2[102400], *put2 = out2;
        ntotalReceived = nlocalsReceived[0];
        local_offsets[0]=nlocalNative; //poure new particles after the OldNatives in local!
        // put2 += sprintf(put2, "%d -> %d part offset %d, ", 0, nlocalsReceived[0], local_offsets[0]);
        for (int p = 1; p<n_proc; p++)
        {
            ntotalReceived += nlocalsReceived[p];
            local_offsets[p] = local_offsets[p-1] + nlocalsReceived[p-1];
            // put2 += sprintf(put2, "%d -> %d part offset %d, ", p, nlocalsReceived[p], local_offsets[p]);
        }
        // printf("I am proc %d and will receive: %s = %d particles\n", rank, out2, ntotalReceived);
        // fflush(stdout);
        MPI_Alltoallv( particles_scatter, partition_sizes, partition_offsets, PARTICLE,
                       local, nlocalsReceived, local_offsets, PARTICLE, MPI_COMM_WORLD );

        if ( (step +1)%SAVEFREQ == 0 and rank == 0)
        {
            time_commu =  read_timer() - time_commu;
            time_binning2 = read_timer();
        }
        //
        //do the binning in each processor of the NEW particles just received.
        //
        for (int i = 0; i<ntotalReceived; i++)  // for the alltoall make sure not to include the native locals!
        {
            // cout << "Proc "<< rank<<" loop Received numb: "<<i<<endl;
            // cout.flush();
            x_bin = floor(local[nlocalNative + i].x/cutoff);
            y_bin = floor(local[nlocalNative + i].y/cutoff);
            bins[ x_bin + binNumPerEdge * y_bin].particlesInBin.push_back( local + nlocalNative + i );
            //update the localNative array by appending the new particles that belong to the processor
            x_proc = x_bin/regColsOfBinsPerProc;
            if (x_proc > colsOfProc -1) x_proc = colsOfProc -1;
            y_proc = y_bin/regRowsOfBinsPerProc;
            if (y_proc > rowsOfProc -1) y_proc = rowsOfProc -1;

            if ( x_proc*rowsOfProc + y_proc == rank )
            {
                localNative[nlocalNativeNew] = nlocalNative + i;
                nlocalNativeNew += 1;
            }
        }     
        if ( (step +1)%SAVEFREQ == 0 and rank == 0)
        {
            time_binning = time_binning + read_timer() - time_binning2;   
            time_computing = read_timer();
        }

        // 
        // print the received and bining of each processor
        // 
        // char out1[102400], *put1 = out1; //fucking learn to clear arrays!!!
        // int partsInBin;
        // for (int b = 0; b < binNumber; b++)
        // {
        //     partsInBin = bins[b].particlesInBin.size();
        //     if (partsInBin > 0)
        //     {
        //         put1 += sprintf(put1, "\tBin %d has %d parts: ", b, partsInBin);
        //         for(int p =0; p< partsInBin; p++)
        //         {
        //             current_part = *(bins[b].particlesInBin[p]);
        //             put1 += sprintf(put1, "%f, ", current_part.x);
        //         }
        //     }
        // }
        // printf("\nI am proc %d and I have received %d part. Now I have %d natives, all bined as:\n %s \n",rank,ntotalReceived, nlocalNativeNew, out1);
        // fflush(stdout);

        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        // if( find_option( argc, argv, "-no" ) == -1 )
        //   if( fsave && (step%SAVEFREQ) == 0 )
        //     save( fsave, nlocalNativeNew, local );

        //
        // Let's apply the forces!
        //
        for (int i =0; i < localNumOfBins; i++)
        {
            apply_force_per_bin(bins, native_bins[i], binNumPerEdge, &dmin, &davg, &navg);
        }
        if ( (step +1)%SAVEFREQ == 0 and rank == 0)
        {
            time_computing = read_timer() - time_computing;
            time_moving = read_timer();
        }

        // printf("Porc %d has FORCE his particles\n",rank);
        // fflush(stdout);
        if( find_option( argc, argv, "-no" ) == -1 )
        {
          
          MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);

 
          if (rank == 0)
          {
            //
            // Computing statistical data
            //
            if (rnavg) {
              absavg +=  rdavg/rnavg;
              nabsavg++;
            }
            if (rdmin < absmin) absmin = rdmin;
          }
        }


        //
        // Now MOVE and redefine the Native particles, putting them in local\
        //
        for( int i = 0; i < nlocalNativeNew; i++ )
        {
            move( local[ localNative[i] ] );
            current_part = local[ localNative[i] ];
            // printf("Proc %d moving part %d with x=%f\n", rank, i , current_part.x);
            // fflush(stdout);
            local[i] = current_part;
        }
        // 
        //set the new boundary of the native particles (including the ones that after moving are not native anymore).
        // 
        nlocalNative = nlocalNativeNew;

        if ( (step +1)%SAVEFREQ == 0 and rank == 0)
        {
            time_moving = read_timer() - time_moving;
            fprintf(fsave, "%f,%f,%f,%f\n", time_binning, time_commu, time_computing, time_moving);
        }
        // printf("Porc %d has MOVE his particles\n",rank);
        // fflush(stdout);

    }

    simulation_time = read_timer( ) - simulation_time;
    // cout << "*****************rank "<<rank<<" finished"<<"**************"<<endl;
    // cout.flush();

    if (rank == 0) 
    {  
      printf( "\n************\nN = %d, n = %d, simulation time = %g seconds\n**********\n",
                                 n_proc,    n,               simulation_time);

      if( find_option( argc, argv, "-no" ) == -1 )
      {
        if (nabsavg) absavg /= nabsavg;
      
        // -The minimum distance absmin between 2 particles during the run of the simulation
        // -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
        // -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
      
        // -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
      
        printf("absmin = %lf, absavg = %lf", absmin, absavg );
        if (absmin < 0.4) fprintf (fsave,"\n\n\nThe minimum distance is below 0.4 meaning that some particle is not interacting\n\n********\n\n");
        if (absavg < 0.8) fprintf (fsave,"\nThe average distance is below 0.8 meaning that most particles are not interacting");
      }
      fprintf(fsave,"\n");     
    // Printing summary data
       
      if( fsum )
      {
         fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
      }
    }
  
    //
    //  release resources
    //
    if ( fsum )
        fclose( fsum );

    free( local );
    free( particles_scatter );
    free( partition_offsets );
    free( partition_sizes );
    // free( localNative );
    // free( native_bins );

    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
