#include <vector>
#include "world_omp.h"
#include <omp.h>
#include <stdlib.h>

using namespace std;


direction Directions[] = { {-1,-1}, {-1,0}, {-1,1}, {0,-1}, {0,0}, {0,1}, {1,-1}, {1,0}, {1,1} };


world_omp_t::world_omp_t(grid_t outGrid, int numParticles, int numBins)
{
	grid = outGrid;
	particles.resize(numParticles);
	divisions.resize(numBins);
}


void world_omp_t::resetDivision()
{
	relocateAllParticles();

	sort(particles.begin(), particles.end(), ascending);

	int divisionLine = 0;
	for (int i = 0; i < particles.size(); i++)
	{
		while(particles[i].binNumber > divisionLine)
		{
			divisions[divisionLine] = i;
			divisionLine++;
		}
	}

	for (int j = divisionLine; j < divisions.size(); j++)
	{
		divisions[j] = particles.size();
	}
}

void world_omp_t::relocateAllParticles()
{
	#pragma omp parallel for schedule(static)
	for(int i = 0; i < particles.size(); i++)
	{
		particles[i].binNumber = relocateOneParticle(particles[i]);
	}
}

int world_omp_t::relocateOneParticle(particle_t &particle)
{
	int bin_x = particle.x / grid.getBinSize();
	int bin_y = particle.y / grid.getBinSize();
	int ans = bin_x + bin_y * grid.getNumBinPerEdge();
	return ans;
}

void world_omp_t::computeAllParticleForces(double *dmin, double *davg, int *navg)
{
       	int numThreads = omp_get_num_threads();
	int force_batch_size = (particles.size() - 1 + numThreads) / numThreads;
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < particles.size(); i+= force_batch_size)
	{
		double individual_dmin = 1.0;
		double individual_davg = 0.0;
		int individual_navg = 0;
		for (int j = i; j < min(i + force_batch_size, particles.size()); j++)
		{
			computeOneParticleForces(particles[j], &individual_dmin, &individual_davg, &individual_navg);
		}

		#pragma omp critical
		{
			*dmin = min(*dmin, individual_dmin);
			*davg += individual_davg;
			*navg += individual_navg;
		}
	}
	
}

void world_omp_t::computeOneParticleForces(particle_t &particle, double *dmin, double *davg, int *navg)
{
	int numBinPerEdge = grid.getNumBinPerEdge();

	int bin_x = particle.binNumber % numBinPerEdge;
	int bin_y = particle.binNumber / numBinPerEdge;

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
			computeOneParticleForcesByBin(particle, binIndex, dmin, davg, navg);
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
				computeOneParticleForcesByBin(particle, binIndex, dmin, davg, navg);
			}
		}
	}
}

void world_omp_t::computeOneParticleForcesByBin(particle_t &particle, int binIndex, double *dmin, double *davg, int *navg)
{
	int start, end;
	//The divisions is like {3, 5, 8, 10 ... particles.size()}
	if (binIndex == 0)
	{
		start = 0;
		end = divisions[0];
	}
	else
	{
		start = divisions[binIndex - 1];
		end = divisions[binIndex];
	}

	for (int i = start; i < end; i++)
	{
		if(!ifSameParticle(particle, particles[i]))
		{
			apply_force(particle, particles[i], dmin, davg, navg);
		}
	}

}

void world_omp_t::moveParticles()
{	
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
	return (a.binNumber == b.binNumber) && (a.x == b.x) && (a.y == b.y);  
}

bool ascending(const particle_t &a, const particle_t &b)
{
	return (a.binNumber < b.binNumber);
}
