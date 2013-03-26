/* ---------------------------------------------------
   FILE:     simParams.h
	AUTHOR:   Josh Bongard
	DATE:     October 5, 2000
	FUNCTION: This class contains all miscellaneous
				 data and functions for this simulation.
 -------------------------------------------------- */

#include "iostream"
#include "fstream"
#include <time.h>

#ifndef _SIM_PARAMS_H
#define _SIM_PARAMS_H

#include "matrix.h"

class SIM_PARAMS {

public:
	int		randSeed;
	double		mutationRate;
	int		numHiddenNeurons;

	char		bodyParamsFileName[200];
	char		objsConqFileName[200];

	time_t		startTime;
	clock_t		startingClock;
	int		evaluationLength;

	int		writeAllReplacements;
	int		evolveLengths;
	int		evolveRadii;
	int		evolveSpacings;

	int		selectForCoreDists;
	int		selectForLeftShoulderDists;
	int		selectForRightShoulderDists;
	MATRIX		*objectOrder;
	int		loadState;

	long		totalTimeStepsUsed;
	int		popSize;
	int		numSlots;
	int		dataFileBufferWidth;
	int		useEarlyStopping;
	int		numLayers;
	int		numIndsPerLayer;
	int		successOccurred;
	int		successfulGenomeIndex;
	int		compress;
	int		noScaffolding;
	int		betweenScaffolding;
	int		withinScaffolding;
	int		withinScaffolding2;
	int		withinScaffolding3;
	int		withinScaffolding4;
	int		totalObjects;
	int		targetDistance;
	int		runIndex;
	int		stanceIndex;
	double		turningAmount;

public:
	SIM_PARAMS(int argc, char **argv);
	~SIM_PARAMS(void);
	void   CloseDataFiles(void);
	void   DirectoryMake(char *dirName);
	void   FileCreate(char *fileName);
	void   FileDelete(char *fileName);
	int    FileExists(char *fileName);
	void   FileRename(char *src, char *dest);
	int    FlipCoin(void);
	ofstream *GetOutFile(char *fileName);
	double HoursSinceStart(void);
	double Min(double x, double y, double z);
	void   ParseParameters(int argc, char **argv);
	void   PermuteList(int *list, int length);
	void   Pause(double s);
	double Rand(double min, double max);
	double RandGaussian(void);
	int    RandInt(int min, int max);
	double Scale(double value, double min1, double max1,
		         double min2, double max2);
	void   SendKillFile(void);
	int    TimeElapsed(void);
	void   WaitForFile(char *fileName);
};

#endif
