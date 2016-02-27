#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "libomp.h"
#include <vector>

using namespace std;

direction Directions[] = {{-1,-1},{-1,0},{-1,1},{0,-1},{0,0},{0,1},{1,-1},{1,0},{1,1}};

void apply_force_per_bin(bin_t* bins, int binIdx, int binNumPerEdge, double* dmin, double* davg,int* navg)
{
    double* dmin1 = dmin;
    double* davg1 = davg;
    int* navg1 = navg;
    for(int i = 0; i < bins[binIdx].particlesInBin.size(); i++)
        apply_force_per_particle(*(bins[binIdx].particlesInBin[i]), bins, binIdx, binNumPerEdge, dmin1, davg1, navg1);
}

void apply_force_per_particle(particle_t &particle, bin_t* bins, int binIdx, int binNumPerEdge, double* dmin, double* davg, int* navg)
{
    particle.ax = particle.ay = 0;
    double* dmin1 = dmin;
    double* davg1 = davg;
    int* navg1 = navg;

    int bin_x = binIdx % binNumPerEdge;
    int bin_y = binIdx / binNumPerEdge;
    
    for(int i = 0; i < 9; i++)
    {
        int bin_x_new = bin_x + Directions[i].x_change;
        int bin_y_new = bin_y + Directions[i].y_change;
        if(bin_x_new >= 0 && bin_x_new < binNumPerEdge && bin_y_new >= 0 && bin_y_new < binNumPerEdge)
        {
            int binIdx_new = bin_x_new + bin_y_new * binNumPerEdge;
            for(int j = 0; j < bins[binIdx_new].particlesInBin.size(); j++)
            {
                apply_force(particle, *(bins[binIdx_new].particlesInBin[j]), dmin1, davg1, navg1);
            }
        }
    }
}
