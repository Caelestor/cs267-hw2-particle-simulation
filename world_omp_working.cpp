#include <vector>
#include "world_omp.h"
#include <stdlib.h>

using namespace std;

#define RELOCATE_BATCH_SIZE	125
#define FORCE_BATCH_SIZE	250
#define MOVE_BATCH_SIZE		125


direction Directions[] = { {-1,-1}, {-1,0}, {-1,1}, {0,-1}, {0,0}, {0,1}, {1,-1}, {1,0}, {1,1} };


world_omp_t::world_omp_t(grid_t outGrid, int numParticles, int numBins)
{
	grid = outGrid;
	particles.resize(numParticles);
	divisions.resize(numBins);
}


void world_omp_t::resetDivision(vector<bin_t> &bins, omp_lock_t* locks)
{
	
	#pragma omp parallel for schedule(static)
	for(int i = 0; i < bins.size(); i++)
		bins[i].particlesInBin.clear();
	
	relocateAllParticles(bins, locks);
}

void world_omp_t::relocateAllParticles(vector<bin_t> &bins, omp_lock_t* locks)
{	
	#pragma omp parallel for schedule(static)
	for(int i = 0; i < particles.size(); i++)
	{
		relocateOneParticle(bins, particles[i], locks);
	}
}

int world_omp_t::relocateOneParticle(vector<bin_t> &bins, particle_t &particle, omp_lock_t* locks)
{
	int bin_x = particle.x / grid.getBinSize();
	int bin_y = particle.y / grid.getBinSize();
	int binNumber = bin_x + bin_y * grid.getNumBinPerEdge();
	particle_t* p = &particle;
	omp_set_lock(locks+binNumber);
	bins[binNumber].particlesInBin.push_back(p);
	omp_unset_lock(locks+binNumber);
}

void world_omp_t::computeAllParticleForces(double *dmin, double *davg, int *navg, vector<bin_t> &bins)
{

	int numThreads = omp_get_num_threads();
	int force_batch_size = (particles.size() - 1 + numThreads) / numThreads;
	#pragma omp parallel for schedule(static)
        for (int i = 0; i < particles.size(); i+=force_batch_size)
	{
		double individual_dmin = 1.0;
		double individual_davg = 0.0;
		int individual_navg = 0;
		vector<bin_t> tmp = bins;
		for (int j = i; j < min(i+force_batch_size, particles.size()); j++)
		{
			computeOneParticleForces(particles[j], &individual_dmin, &individual_davg, &individual_navg, tmp);
		}
 	/*	#pragma omp critical
		{
			*dmin = min(*dmin, individual_dmin);
			*davg += individual_davg;
			*navg += individual_navg;
		}*/
	}
	/*for (int i = 0; i < particles.size(); i += 1)
	{
		double tmp0 = 1.0;
		double tmp1 = 0.0;
		int tmp2 = 0;
		vector<bin_t> tmp = bins;
		computeOneParticleForces(particles[i], &tmp0, &tmp1, &tmp2, tmp);
	}*/
}

void world_omp_t::computeOneParticleForces(particle_t &particle, double *dmin, double *davg, int *navg, vector<bin_t> &bins)
{
	int numBinPerEdge = grid.getNumBinPerEdge();

	int bin_x = particle.x / grid.getBinSize();
	int bin_y = particle.y / grid.getBinSize();

	particle.ax = 0;
	particle.ay = 0;

	//When the particle belongs to a bin that is not an edge bin, i.e., the bin has eight neighbor bins:
	if(numBinPerEdge > 2 && (bin_x > 0 && bin_x < numBinPerEdge - 1) && (bin_y > 1 && bin_y < numBinPerEdge - 1))
	{
		for(int i = 0; i < 9; i += 1)
		{
			int bin_x1 = bin_x + Directions[i].x_change;
			int bin_y1 = bin_y + Directions[i].y_change;
			int binIndex = bin_x1 + bin_y1 * numBinPerEdge;
			computeOneParticleForcesByBin(particle, binIndex, dmin, davg, navg, bins);
		}
	}
	//When the particle belongs to a bin that is an edge bin, i.e., the bin has less than eight neighbor bins:
	else
	{
		for(int i = 0; i < 9; i += 1)
		{
			int bin_x1 = bin_x + Directions[i].x_change;
			int bin_y1 = bin_y + Directions[i].y_change;
			if((bin_x1 >= 0 && bin_x1 < numBinPerEdge) && (bin_y1 >= 0 && bin_y1 < numBinPerEdge))
			{
				int binIndex = bin_x1 + bin_y1 * numBinPerEdge;
				computeOneParticleForcesByBin(particle, binIndex, dmin, davg, navg, bins);
			}
		}
	}
}

void world_omp_t::computeOneParticleForcesByBin(particle_t &particle, int binIndex, double *dmin, double *davg, int *navg, vector<bin_t> &bins)
{
	int iterNum = bins[binIndex].particlesInBin.size();
	for(int i = 0; i < iterNum; i++)
	{
		apply_force(particle, *(bins[binIndex].particlesInBin[i]), dmin, davg, navg);
	}

}

void world_omp_t::moveParticles()
{	
	//Parallelizable
	#pragma omp parallel for schedule(static)
	for(int i = 0; i < particles.size(); i++)
	{
		move(particles[i]);
	}
}


/*vector<int> world_omp_t::getDivisions()
{
	return divisions;
}*/

grid_t world_omp_t::getGrid()
{
	return grid;
}

void world_omp_t::setDivisions(vector<int> new_divisions)
{
	divisions = new_divisions;
}

void world_omp_t::setGrid(grid_t new_grid)
{
	grid = new_grid;
}

bool ifSameParticle(const particle_t &a, const particle_t &b)
{
	return (a.x == b.x) && (a.y == b.y);  
}
