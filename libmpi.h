#ifndef LIBOMP_H
#define LIBOMP_H

#include <vector>
#include <omp.h>
#include "common_mpi.h"
#include <stdlib.h>

using namespace std;

class bin_t
{
public: 
    vector<particle_t*> particlesInBin;

    bin_t(){particlesInBin.resize(100);}

};

struct direction
{
	int x_change;
	int y_change;
};

void apply_force_per_bin(bin_t* bins, int binIdx, int binNumPerEdge, double* dmin, double *davg, int *navg);

void apply_force_per_particle(particle_t &particle, bin_t* bins, int binIdx, int binNumPerEdge, double* dmin, double *davg, int *navg);





#endif
