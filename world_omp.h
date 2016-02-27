#ifndef WORLD_H
#define WORLD_H

#include <vector>
#include "grid.h"
#include <stdlib.h>
#include <algorithm>
#include <omp.h>

using namespace std;

struct direction
{
	int x_change;
	int y_change;
};

class world_omp_t {

private:
	grid_t grid;//The topology of the world
	//The particles inside the world
	vector<int> divisions;//The division lines of the particles. For example: divisions = {0,0,3,5,5,5,10,12,13,16,16} means
	//divisions[0] = 0: [particle 0,particle 0) 	belongs to division 0 (no particles belong to division 0);
	//divisions[1] = 0: [particle 0,particle 0) 	belongs to division 1 (no particles belong to division 1);
	//divisions[2] = 3: [particle 0,particle 3) 	belongs to division 2;
	//divisions[3] = 5: [particle 3,particle 5) 	belongs to division 3;
	//divisions[4] = 5: [particle 5,particle 5) 	belongs to division 4 (no particles belong to division 4);
	//divisions[5] = 5: [particle 5,particle 5) 	belongs to division 5 (no particles belong to division 5);
	//divisions[6] = 10: [particle 5,particle 10) 	belongs to division 6;
	//divisions[7] = 12: [particle 10,particle 12) 	belongs to division 7;
	//divisions[8] = 13: [particle 12,particle 13) 	belongs to division 8;
	//divisions[9] = 16: [particle 13,particle 16) 	belongs to division 9;
	//divisions[10] = 16: [particle 16,particle 16) belongs to division 10 (no particles belong to division 10);
	//Therefore, particles = {2,2,2,3,3,6,6,6,6,6,7,7,8,9,9,9}

public:
  vector<particle_t> particles;
  
	world_omp_t(grid_t outGrid, int numParticles, int numBins);

	void resetDivision();

	void relocateAllParticles();

	int relocateOneParticle(particle_t &p);

	void computeAllParticleForces(double *dmin, double *davg, int *navg);

	void computeOneParticleForces(particle_t &particle, double *dmin, double *davg, int *navg);

	void computeOneParticleForcesByBin(particle_t &particle, int binIndex, double *dmin, double *davg, int *navg);

	void moveParticles();

//	vector<int> getDivisions();

	grid_t getGrid();

	void setDivisions(vector<int> new_divisions);

	void setGrid(grid_t new_grid);

};

bool ifSameParticle(const particle_t &a, const particle_t &b);

bool ascending(const particle_t &a, const particle_t &b);

#endif
