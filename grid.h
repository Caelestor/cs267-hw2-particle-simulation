#ifndef GRID_H
#define GRID_H

#include <vector>
#include <cstdio>
#include "common_new.h"


class grid_t {
private:
	double gridSize;
	double binSize;
	int numParticle;
	int numBinPerEdge;

public:
  grid_t();
	grid_t(double myGridSize, double myBinSize, int myParticleNum);

	double getGridSize();
	double getBinSize();
	int getNumParticle();
	int getNumBinPerEdge();
	int getNumBin();

};

#endif


