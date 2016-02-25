#include <vector>
#include "world.h"
#include <stdlib.h>

using namespace std;


direction Directions[] = { {-1,-1}, {-1,0}, {-1,1}, {0,-1}, {0,0}, {0,1}, {1,-1}, {1,0}, {1,1} };


world_t::world_t(grid_t outGrid, int numParticles, int numBins)
{
	grid = outGrid;
	particles.resize(numParticles);
	divisions.resize(numBins);
}


void world_t::resetDivision()
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

void world_t::relocateAllParticles()
{
	//Parallalizable
	for(auto & particle : particles)
	{
		particle.binNumber = relocateOneParticle(particle);
	}
}

int world_t::relocateOneParticle(particle_t &particle)
{
	int bin_x = particle.x / grid.getBinSize();
	int bin_y = particle.y / grid.getBinSize();
	int ans = bin_x + bin_y * grid.getNumBinPerEdge();
	return ans;
}

void world_t::computeAllParticleForces(double *dmin, double *davg, int *navg)
{
	for (auto & particle : particles)
	{
		computeOneParticleForces(particle, dmin, davg, navg);
	}
}

void world_t::computeOneParticleForces(particle_t &particle, double *dmin, double *davg, int *navg)
{
	int numBinPerEdge = grid.getNumBinPerEdge();

	int bin_x = particle.binNumber % numBinPerEdge;
	int bin_y = particle.binNumber / numBinPerEdge;

	particle.ax = 0;
	particle.ay = 0;

	//When the particle belongs to a bin that is not an edge bin, i.e., the bin has eight neighbor bins:
	if(numBinPerEdge > 2 && (bin_x >= 1 && bin_x <= numBinPerEdge - 1) && (bin_y >= 1 && bin_y <= numBinPerEdge - 1))
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
			if((bin_x1 >= 0 && bin_x1 <= numBinPerEdge) && (bin_y1 >= 0 && bin_y1 <= numBinPerEdge))
			{
				int binIndex = bin_x1 + bin_y1 * numBinPerEdge;
				computeOneParticleForcesByBin(particle, binIndex, dmin, davg, navg);
			}
		}
	}
}

void world_t::computeOneParticleForcesByBin(particle_t &particle, int binIndex, double *dmin, double *davg, int *navg)
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

void world_t::moveParticles()
{	
	//Parallelizable
	for(auto & particle: particles)
	{
		move(particle);
	}
}


/*vector<int> world_t::getDivisions()
{
	return divisions;
}*/

grid_t world_t::getGrid()
{
	return grid;
}

void world_t::setDivisions(vector<int> new_divisions)
{
	divisions = new_divisions;
}

void world_t::setGrid(grid_t new_grid)
{
	grid = new_grid;
}

bool ifSameParticle(particle_t &a, particle_t &b)
{
	return (a.binNumber == b.binNumber) && (a.x == b.x) && (a.y == b.y);  
}

bool ascending(particle_t &a, particle_t &b)
{
	return (a.binNumber < b.binNumber);
}
