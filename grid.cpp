#include <vector>
#include <cmath>
#include <cstdio>
#include "grid.h"


grid_t::grid_t(double myGridSize, double myBinSize, int myParticleNum)
{
	gridSize = myGridSize;
	binSize = myBinSize;
	numParticle = myParticleNum;
	if(myGridSize < myBinSize)
	{
		binSize = gridSize;
		numBinPerEdge = 1;
	}
	else
		numBinPerEdge = ceil(myGridSize / myBinSize);
	//Do we need to consider if the gridSize is not a multiple of binSize?

}

grid_t::grid_t(){}

double grid_t::getGridSize()
{
	return gridSize;
}

double grid_t::getBinSize()
{
	return binSize;
}

int grid_t::getNumParticle()
{
	return numParticle;
}

int grid_t::getNumBinPerEdge()
{
	return numBinPerEdge;
}

int grid_t::getNumBin()
{
	return numBinPerEdge * numBinPerEdge;
}

