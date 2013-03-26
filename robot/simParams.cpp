/* ---------------------------------------------------
   FILE:     simParams.cpp
	AUTHOR:   Josh Bongard
	DATE:     October 5, 2000
	FUNCTION: This class contains all miscellaneous
				 data and functions for this simulation.
 -------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "sys/stat.h"

#ifndef _SIM_PARAMS_CPP
#define _SIM_PARAMS_CPP

#include "simParams.h"

extern int		RANDOM_SEED;
extern char		DATA_DIRECTORY[100];
extern int		DATA_FILE_BUFFER;
extern int		NUM_MUTATIONS;
extern int		NUM_HIDDEN_NEURONS;

extern char		TEMP_FILENAME_PREFIX_NEW[100];
extern char		TEMP_FILENAME[100];
extern char		BODY_OUT_FILENAME[100];
extern char		SENSOR_IN_FILENAME[100];

extern int		EVALUATION_LENGTH;
extern int		TARGET_DISTANCE;

extern double		DEFAULT_MUTATION_RATE;
extern int		NUM_SLOTS;
extern int		TOTAL_OBJECTS;

SIM_PARAMS::SIM_PARAMS(int argc, char **argv) {

	startTime = time(NULL);
	startingClock = clock();

	randSeed = RANDOM_SEED;
	mutationRate = DEFAULT_MUTATION_RATE;
	numHiddenNeurons = NUM_HIDDEN_NEURONS;

	loadState = false;

	evaluationLength    = EVALUATION_LENGTH;
	targetDistance	    = TARGET_DISTANCE;
	turningAmount	    = 0.0;

	noScaffolding	   = true;
	betweenScaffolding = false;
	withinScaffolding = false;
	withinScaffolding2 = false;
	withinScaffolding3 = false;
	withinScaffolding4 = false;

	totalTimeStepsUsed  = 0;

	writeAllReplacements = false;
	evolveLengths = false;
	evolveRadii = false;
	evolveSpacings = false;

	selectForCoreDists		= false;
	selectForLeftShoulderDists			= false;
	selectForRightShoulderDists		= false;

	numSlots = NUM_SLOTS;
	useEarlyStopping = false;
	numLayers = 1;
	compress = false;

	totalObjects = TOTAL_OBJECTS;
	stanceIndex=0;

	ParseParameters(argc,argv);

	//srand(randSeed);

	numIndsPerLayer = int(floor(400.0/double(numLayers)));

	popSize = numLayers * numIndsPerLayer;

	objectOrder = new MATRIX(1,totalObjects,0.0);
	for (int j=0;j<totalObjects;j++)
		objectOrder->Set(0,j,j);

	/*
	// Randomize object order
	if ( totalObjects>1 )
		for (int j=0;j<1000;j++) {
		        int index1 = RandInt(0,totalObjects-1);
		        int index2 = RandInt(0,totalObjects-1);
		        while( index2 == index1 )
		                index2 = RandInt(0,totalObjects-1);
		        int temp = int(objectOrder->Get(0,index1));
		        objectOrder->Set(0,index1,objectOrder->Get(0,index2));
		        objectOrder->Set(0,index2,temp);
		}
	*/

/*	
	// Select for locomotion, then turning
	int leftPlacement = int(double(totalObjects)/2.0 - 1.0);
	int rightPlacement = int(double(totalObjects)/2.0);
	for (int i=0;i<totalObjects;i=i+2) {
		objectOrder->Set(0,i,leftPlacement);
		objectOrder->Set(0,i+1,rightPlacement);
		leftPlacement--;
		rightPlacement++;
	}
*/
	objectOrder->Print();

	dataFileBufferWidth = 7 + 4*popSize;

	successOccurred = false;
	successfulGenomeIndex = -1;

	runIndex = randSeed - 100*int(double(randSeed)/100.0);
}

SIM_PARAMS::~SIM_PARAMS(void) {

	if ( objectOrder ) {
		delete objectOrder;
		objectOrder = NULL;
	}

	SendKillFile();
}

void  SIM_PARAMS::CloseDataFiles(void) {

}

void  SIM_PARAMS::DirectoryMake(char *dirName) {

	char command[500];

	sprintf(command,"rd %s",dirName);
	system(command);

	sprintf(command,"md %s",dirName);
	system(command);
}

void  SIM_PARAMS::FileCreate(char *fileName) {

}

void  SIM_PARAMS::FileDelete(char *fileName) {

	char command[100];

	while ( FileExists(fileName) ) {
		sprintf(command,"del %s",fileName);
		system(command);
	}
}

int   SIM_PARAMS::FileExists(char *fileName) {

  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  intStat = stat(fileName,&stFileInfo);
  if(intStat == 0) {
    blnReturn = true;
  } else {
    blnReturn = false;
  }
  
  return(blnReturn);
}

void SIM_PARAMS::FileRename(char *src, char *dest) {

	char command[200];

	sprintf(command,"ren %s %s",src,dest);

	system(command);
}

int SIM_PARAMS::FlipCoin(void) {

	return( Rand(0.0,1.0) < 0.5 );
}

ofstream *SIM_PARAMS::GetOutFile(char *fileName) {

	ofstream *outFile = new ofstream(fileName);

	return( outFile );
}

double SIM_PARAMS::HoursSinceStart(void) {

	double hoursSinceStart = double(time(NULL) - startTime) / 3600.0;

	return( hoursSinceStart );
}

double SIM_PARAMS::Min(double x, double y, double z) {

	if ( (x<=y) && (x<=z) )
		return(x);
	else
		if ( (y<=x) && (y<=z) )
			return(y);
		else
			return(z);
}

void  SIM_PARAMS::ParseParameters(int argc, char **argv) {

	int currParam;

	for(currParam=0;currParam<argc;currParam++) {

		if ( strcmp(argv[currParam],"-r") == 0 )
			randSeed = atoi(argv[currParam+1]);

		if ( strcmp(argv[currParam],"-m") == 0 )
			mutationRate = atof(argv[currParam+1]);

		if ( strcmp(argv[currParam],"-h") == 0 )
			numHiddenNeurons = atoi(argv[currParam+1]);

		if ( strcmp(argv[currParam],"-el") == 0 )
			evaluationLength = atoi(argv[currParam+1]);

                if ( strcmp(argv[currParam],"-td") == 0 )
                        targetDistance = atoi(argv[currParam+1]);

                if ( strcmp(argv[currParam],"-ta") == 0 )
                        turningAmount = atof(argv[currParam+1]);

		if ( strcmp(argv[currParam],"-core") == 0 )
			selectForCoreDists = atoi(argv[currParam+1]);

		if ( strcmp(argv[currParam],"-left") == 0 )
			selectForLeftShoulderDists = atoi(argv[currParam+1]);

		if ( strcmp(argv[currParam],"-right") == 0 )
			selectForRightShoulderDists = atoi(argv[currParam+1]);

		if ( strcmp(argv[currParam],"-s") == 0 )
			numSlots = atoi(argv[currParam+1]);

		if ( strcmp(argv[currParam],"-p") == 0 )
			popSize = atoi(argv[currParam+1]);

		if ( strcmp(argv[currParam],"-war") == 0 )
			writeAllReplacements = true;

		if ( strcmp(argv[currParam],"-es") == 0 )
			useEarlyStopping = true;

                if ( strcmp(argv[currParam],"-ua") == 0 )
                        numLayers = atoi(argv[currParam+1]);

		if ( strcmp(argv[currParam],"-co") == 0 )
			compress = true;

		if ( strcmp(argv[currParam],"-load") == 0 )
			loadState = true;

		if ( strcmp(argv[currParam],"-bs") == 0 ) {
			noScaffolding = false;
			betweenScaffolding = true;
		}
		if ( strcmp(argv[currParam],"-ws") == 0 ) {
                        noScaffolding = false;
			withinScaffolding = true;
		}
		if ( strcmp(argv[currParam],"-ws2") == 0 ) {
                        noScaffolding = false;
			withinScaffolding2 = true;
		}
		if ( strcmp(argv[currParam],"-ws3") == 0 ) {
                        noScaffolding = false;
			withinScaffolding3 = true;
		}
		if ( strcmp(argv[currParam],"-ws4") == 0 ) {
                        noScaffolding = false;
			withinScaffolding4 = true;
		}
		if ( strcmp(argv[currParam],"-to") == 0 ) {
			totalObjects = atoi(argv[currParam+1]);
		}
	}
}

void SIM_PARAMS::PermuteList(int *list, int length) {

	int firstIndex;
	int secondIndex;
	int temp;

	for (int i=0;i<1000;i++) {

		firstIndex = RandInt(0,length-1);
		secondIndex = RandInt(0,length-1);

		while ( secondIndex==firstIndex )
			secondIndex = RandInt(0,length-1);

		temp = list[firstIndex];
		list[firstIndex] = list[secondIndex];
		list[secondIndex] = temp;
	}
}

void   SIM_PARAMS::Pause(double s) {

	double pauseStart = double(time(NULL));
	double currTime   = double(time(NULL));
	double timeUntilDone;

	while ( currTime - pauseStart < s ) {
		timeUntilDone = s - (currTime - pauseStart);
		currTime   = double(time(NULL));
	}
	printf("Pause finished.\n");
}

double SIM_PARAMS::Rand(double min, double max) {

	double zeroToOne = ((double)rand()) / RAND_MAX;
	double returnVal;

	returnVal = (zeroToOne * (max-min)) + min;
	return returnVal;
}

double SIM_PARAMS::RandGaussian(void) {

	double w = 1.01;
	double x1, x2;

	while ( w >= 1.0 ) {
		x1 = 2.0*Rand(0,1) - 1;
		x2 = 2.0*Rand(0,1) - 1;
		w = x1*x1 + x2*x2;
	}	
	w = sqrt( (-2.0*log(w))/w );
	
	return( x1*w );
}

int SIM_PARAMS::RandInt(int min, int max) {

	if ( min == max )
		return( min );
	else {
		int val = (rand() % (max-min+1)) + min;
		if ( val > max )
			val = max;
		if ( val < min )
			val = min;	
		return( val );
	}
}

double SIM_PARAMS::Scale(double value, double min1, double max1,
								 double min2, double max2) {

	if ( min1 < 0 )
		value = value - min1;
	else
		value = value + min1;

	return( (value*(max2-min2)/(max1-min1)) + min2 );
}

void SIM_PARAMS::SendKillFile(void) {

	char killFileName[200];

	for (int i=0;i<numSlots;i++) {

		sprintf(killFileName,	"%s%d_%d/Kill.dat",TEMP_FILENAME_PREFIX_NEW,randSeed,i);
		ofstream *killFile = new ofstream(killFileName);
		killFile->close();
		delete killFile;
		killFile = NULL;
	}
}

int  SIM_PARAMS::TimeElapsed(void) {

/*
	double minutesElapsed = HoursSinceStart()*60.0;
	double secondsElapsed = minutesElapsed*60.0;

	return( secondsElapsed > 10.0 );
*/

	return( false );
}

void SIM_PARAMS::WaitForFile(char *fileName) {

	while ( !FileExists(fileName) );
}

#endif
