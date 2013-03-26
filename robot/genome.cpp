/* ---------------------------------------------------
   FILE:     genome.cpp
	AUTHOR:   Josh Bongard
	DATE:     October 2, 2000
	FUNCTION: This class contains all information for
				 a single physical segment of an organism
				 in the MathEngine environment.
 -------------------------------------------------- */

#include <math.h>

#ifndef _GENOME_CPP
#define _GENOME_CPP

#include "genome.h"
#include "matrix.h"
#include "simParams.h"

extern SIM_PARAMS *simParams;
extern char SENSOR_IN_FILENAME[100];

extern double TAU_MIN;
extern double TAU_MAX;
extern double WEIGHT_MIN;
extern double WEIGHT_MAX;
extern double OMEGA_MIN;
extern double OMEGA_MAX;
extern double SENSOR_MIN;
extern double SENSOR_MAX;

extern int    NUM_GENERATIONS;

extern int    NUM_RAY_SENSORS;
extern int    NUM_TOUCH_SENSORS;
extern int    NUM_ANGLE_SENSORS;
extern int    NUM_SENSORS;
extern int    NUM_NEURONS;

extern double    NUM_FINGERS;
extern double    FINGER_LENGTH;
extern double	 MIN_FINGER_LENGTH;
extern double	 MAX_FINGER_LENGTH;
extern double	 FINGER_RADIUS;
extern double	 MIN_FINGER_RADIUS;
extern double	 MAX_FINGER_RADIUS;
extern double	 TARGET_RADIUS;

extern double	 CORE_DIST_MAX;
extern double    LEFT_SHOULDER_DIST_MAX;
extern double    RIGHT_SHOULDER_DIST_MAX;

extern char	TAUS_OUT_FILENAME[100];
extern char	WEIGHTS_OUT_FILENAME[100];
extern char	OMEGAS_OUT_FILENAME[100];
extern char	SENSORS_OUT_FILENAME[100];

GENOME::GENOME(int evalsCurrent) {

	totalNeurons = NUM_NEURONS+simParams->numHiddenNeurons;

	taus = new MATRIX(1,totalNeurons);
	for (int i=0;i<totalNeurons;i++)
		taus->Set(0,i,simParams->Rand(TAU_MIN,TAU_MAX));

	weights = new MATRIX(totalNeurons,totalNeurons);
	for (int i=0;i<totalNeurons;i++)
		for (int j=0;j<totalNeurons;j++)
			weights->Set(i,j,simParams->Rand(WEIGHT_MIN,WEIGHT_MAX));

	omegas = new MATRIX(1,totalNeurons);
	for (int i=0;i<totalNeurons;i++)
		omegas->Set(0,i,simParams->Rand(OMEGA_MIN,OMEGA_MAX));

	sensorWeights = new MATRIX(NUM_SENSORS,totalNeurons);
	for (int i=0;i<NUM_SENSORS;i++)
		for (int j=0;j<totalNeurons;j++)
			sensorWeights->Set(i,j,simParams->Rand(SENSOR_MIN,SENSOR_MAX));

	fingerLengths = new MATRIX(3,int(NUM_FINGERS),FINGER_LENGTH);
	fingerRadii   = new MATRIX(3,int(NUM_FINGERS),FINGER_RADIUS);
	fingerSpacings = new MATRIX(1,int(NUM_FINGERS),0.0);
	for (int j=0;j<NUM_FINGERS;j++)
		fingerSpacings->Set(0,j,j*2.0*3.14159/NUM_FINGERS);

	innerFingerRanges = new MATRIX(2,int(NUM_FINGERS),0.0);
	for (int j=0;j<NUM_FINGERS;j++) {
		innerFingerRanges->Set(0,j,-0.01);
		innerFingerRanges->Set(1,j,+0.01);
	}

	objectsToEvaluate = 1;

	Initialize();

	childDeathDistribution = new MATRIX(1,simParams->totalObjects,0.0);

	myParent = NULL;

	evalsCreated 	= evalsCurrent;
}

GENOME::GENOME(GENOME *parent) {

	totalNeurons = parent->totalNeurons;

	taus = new MATRIX(parent->taus);

	weights = new MATRIX(parent->weights);

	omegas = new MATRIX(parent->omegas);

	sensorWeights = new MATRIX(parent->sensorWeights);

	fingerLengths = new MATRIX(parent->fingerLengths);
	fingerRadii   = new MATRIX(parent->fingerRadii);
	fingerSpacings = new MATRIX(parent->fingerSpacings);
	innerFingerRanges = new MATRIX(parent->innerFingerRanges);

	objectsToEvaluate = parent->objectsToEvaluate;

	Initialize();

	childDeathDistribution = new MATRIX(parent->childDeathDistribution);
	myParent = parent;

	evalsCreated = parent->evalsCreated;
}

GENOME::~GENOME(void) {

	Clear();

	if ( childDeathDistribution ) {
		delete childDeathDistribution;
		childDeathDistribution = NULL;
	}

	delete innerFingerRanges;
	innerFingerRanges = NULL;

	delete fingerSpacings;
	fingerSpacings = NULL;

	delete fingerLengths;
	fingerLengths = NULL;

	delete fingerRadii;
	fingerRadii = NULL;

	delete taus;
	taus = NULL;

	delete weights;
	weights = NULL;

	delete omegas;
	omegas = NULL;

	delete sensorWeights;
	sensorWeights = NULL;
}

int  GENOME::Age(int evalsCurrent) {

	double age = evalsCurrent - evalsCreated;
	age = age / double(simParams->popSize);
	age = 1.0 + age;

	return( int(floor(age)) );
}

int GENOME::AllObjectsEvaluated(void) {

	return( objectsEvaluated == objectsToEvaluate );
}

void GENOME::Balance(int objectsOnLayer) {

	while ( objectsToEvaluate > objectsOnLayer )

		Compress();

	while ( objectsToEvaluate < objectsOnLayer )

		Expand();
}

void GENOME::CalculateFitness(	int currentGeneration, int rowIndex,
				int currentSlot) {

	char sensorFileName[200];
	sprintf(sensorFileName,"/tmp/Files%d_%d/%s",simParams->randSeed,currentSlot,SENSOR_IN_FILENAME);
	ifstream *inFile = new ifstream(sensorFileName,ios::in);

	double tmp;
	(*inFile) >> tmp;

	int currEvalLength;
	(*inFile) >> currEvalLength;

	evalLengths->Set(0,objectsEvaluated,currEvalLength);
	simParams->totalTimeStepsUsed = simParams->totalTimeStepsUsed + int(evalLengths->Get(0,objectsEvaluated));

	inFile->close();
	delete inFile;
	inFile = NULL;
	simParams->FileDelete(sensorFileName);

       if ( (tmp==0) || (currEvalLength==1) ) {

		Invalidate();

		char fileName[200];
		sprintf(fileName,"/tmp/Files%d_%d/myLeftShoulderDists.dat",simParams->randSeed,currentSlot);
		simParams->FileDelete(fileName);

		return;
	}

	CalculateCoreDist(          objectsEvaluated,currentSlot);
	CalculateLeftShoulderDist(  objectsEvaluated,currentSlot);
	CalculateRightShoulderDist( objectsEvaluated,currentSlot);

	objectsEvaluated++;

	fitness = coreDist * leftShoulderDist * rightShoulderDist;
}

void   GENOME::Compress(void) {

	MATRIX *temp;

	objectsToEvaluate--;

	if ( objectsEvaluated > objectsToEvaluate )
		objectsEvaluated--;

	temp = coreDists;
	coreDists = new MATRIX(1,objectsToEvaluate,CORE_DIST_MAX);
	for (int j=0;j<objectsToEvaluate;j++)
		coreDists->Set(0,j,temp->Get(0,j));
	delete temp;
	temp = NULL;
	if ( simParams->selectForCoreDists )
		coreDist = coreDists->Mean();
	else
		coreDist = CORE_DIST_MAX;

	temp = leftShoulderDists;
	leftShoulderDists = new MATRIX(1,objectsToEvaluate,LEFT_SHOULDER_DIST_MAX);
	for (int j=0;j<objectsToEvaluate;j++)
		leftShoulderDists->Set(0,j,temp->Get(0,j));
	delete temp;
	temp = NULL;
	if ( simParams->selectForLeftShoulderDists )
		leftShoulderDist = leftShoulderDists->Mean();
	else
		leftShoulderDist = LEFT_SHOULDER_DIST_MAX;

	temp = rightShoulderDists;
	rightShoulderDists = new MATRIX(1,objectsToEvaluate,RIGHT_SHOULDER_DIST_MAX);
	for (int j=0;j<objectsToEvaluate;j++)
		rightShoulderDists->Set(0,j,temp->Get(0,j));
	delete temp;
	temp = NULL;
	if ( simParams->selectForRightShoulderDists )
		rightShoulderDist = rightShoulderDists->Mean();
	else
		rightShoulderDist = RIGHT_SHOULDER_DIST_MAX;

	delete objectOrdering;
	objectOrdering 		= new MATRIX(1,objectsToEvaluate,0.0);
	for (int j=0;j<objectsToEvaluate;j++)
		objectOrdering->Set(0,j,simParams->objectOrder->Get(0,j));

	temp = evalLengths;
	evalLengths = new MATRIX(1,objectsToEvaluate,0.0);
	for (int j=0;j<objectsToEvaluate;j++)
		evalLengths->Set(0,j,temp->Get(0,j));
	delete temp;
	temp = NULL;

	fitness = coreDist * leftShoulderDist * rightShoulderDist;
}

int GENOME::Dominates(GENOME *otherGenome) {

	if ( otherGenome->Successful() )
		return( false );

	//RISI replaced round with floor
	int dominatesCore 	= 	int(std::floor(1000000.0*coreDist)) >= 
					int(std::floor(1000000.0*otherGenome->coreDist));

	int dominatesLeft 	= 	int(std::floor(1000000.0*leftShoulderDist)) >=
					int(std::floor(1000000.0*otherGenome->leftShoulderDist));

	int dominatesRight 	= 	int(std::floor(1000000.0*rightShoulderDist)) >= 
					int(std::floor(1000000.0*otherGenome->rightShoulderDist));

	return( dominatesCore && dominatesLeft && dominatesRight );
}

void   GENOME::Expand(void) {

	MATRIX *temp;

	objectsToEvaluate++;

	temp = coreDists;
	coreDists = new MATRIX(1,objectsToEvaluate,CORE_DIST_MAX);
	for (int j=0;j<objectsToEvaluate-1;j++)
		coreDists->Set(0,j,temp->Get(0,j));
	delete temp;
	temp = NULL;
	if ( simParams->selectForCoreDists )
		coreDist = coreDists->Mean();
	else
		coreDist = CORE_DIST_MAX;

	temp = leftShoulderDists;
	leftShoulderDists = new MATRIX(1,objectsToEvaluate,LEFT_SHOULDER_DIST_MAX);
	for (int j=0;j<objectsToEvaluate-1;j++)
		leftShoulderDists->Set(0,j,temp->Get(0,j));
	delete temp;
	temp = NULL;
	if ( simParams->selectForLeftShoulderDists )
		leftShoulderDist = leftShoulderDists->Mean();
	else
		leftShoulderDist = LEFT_SHOULDER_DIST_MAX;

	temp = rightShoulderDists;
	rightShoulderDists = new MATRIX(1,objectsToEvaluate,RIGHT_SHOULDER_DIST_MAX);
	for (int j=0;j<objectsToEvaluate-1;j++)
		rightShoulderDists->Set(0,j,temp->Get(0,j));
	delete temp;	
	temp = NULL;
	if ( simParams->selectForRightShoulderDists )
		rightShoulderDist = rightShoulderDists->Mean();
	else
		rightShoulderDist = RIGHT_SHOULDER_DIST_MAX;

	delete objectOrdering;
	objectOrdering 		= new MATRIX(1,objectsToEvaluate,0.0);
	for (int j=0;j<objectsToEvaluate;j++)
		objectOrdering->Set(0,j,simParams->objectOrder->Get(0,j));

	temp = evalLengths;
	evalLengths = new MATRIX(1,objectsToEvaluate,0.0);
	for (int j=0;j<objectsToEvaluate-1;j++)
		evalLengths->Set(0,j,temp->Get(0,j));
	delete temp;
	temp = NULL;

	fitness = coreDist * leftShoulderDist * rightShoulderDist;
}

double GENOME::MeanEvalTimes(void) {

	double mean = 0.0;
	for (int j=0;j<objectsEvaluated;j++)
		mean = mean + evalLengths->Get(0,j);

	if ( objectsEvaluated==0 )
		return( evalLengths->Get(0,0) );
	else
		return( mean/double(objectsEvaluated) );
}

void GENOME::Mutate(void) {

	int numVals = 	totalNeurons + 
			totalNeurons*totalNeurons + 
			totalNeurons + 
			NUM_SENSORS*totalNeurons;

	double mutProb = simParams->mutationRate;
	int mutationMade = false;
	double mutChange;

	while ( !mutationMade ) {

		for (int j=0;j<totalNeurons;j++)
			if ( simParams->Rand(0,1) < mutProb ) {
				mutChange = simParams->RandGaussian(); // between -3 and 3
				mutChange = mutChange + 3.0; // between 0 and 6
				mutChange = mutChange / 6.0; // between 0 and 1
				mutChange = mutChange * (TAU_MAX - TAU_MIN); // between 0 and max-min
				mutChange = mutChange + TAU_MIN; // between min and max 		
				taus->Set(0,j,mutChange);
				mutationMade = true;
			}

		for (int i=0;i<totalNeurons;i++)
			for (int j=0;j<totalNeurons;j++)
				if ( simParams->Rand(0,1) < mutProb ) {
					mutChange = simParams->RandGaussian(); // between -3 and 3
					mutChange = mutChange + 3.0; // between 0 and 6
					mutChange = mutChange / 6.0; // between 0 and 1
					mutChange = mutChange * (WEIGHT_MAX - WEIGHT_MIN); // between 0 and max-min
					mutChange = mutChange + WEIGHT_MIN; // between min and max 
					weights->Set(i,j,mutChange);
					mutationMade = true;
				}

		for (int j=0;j<totalNeurons;j++)
			if ( simParams->Rand(0,1) < mutProb ) {
				mutChange = simParams->RandGaussian(); // between -3 and 3
				mutChange = mutChange + 3.0; // between 0 and 6
				mutChange = mutChange / 6.0; // between 0 and 1
				mutChange = mutChange * (OMEGA_MAX - OMEGA_MIN); // between 0 and max-min
				mutChange = mutChange + OMEGA_MIN; // between min and max 
				omegas->Set(0,j,mutChange);
				mutationMade = true;
			}

		for (int i=0;i<NUM_SENSORS;i++)
			for (int j=0;j<totalNeurons;j++)
				if ( simParams->Rand(0,1) < mutProb ) {
					mutChange = simParams->RandGaussian(); // between -3 and 3
					mutChange = mutChange + 3.0; // between 0 and 6
					mutChange = mutChange / 6.0; // between 0 and 1
					mutChange = mutChange * (SENSOR_MAX - SENSOR_MIN); // between 0 and max-min
					mutChange = mutChange + SENSOR_MIN; // between min and max 
					sensorWeights->Set(i,j,mutChange);
					mutationMade = true;
				}
	}
}

void GENOME::Print(void) {

	printf("%3.3f \t %3.3f \t %3.3f \n",coreDist,leftShoulderDist,rightShoulderDist);
}

void GENOME::Reset(void) {

	Clear();

	Initialize();
}

void GENOME::Save(ofstream *outFile) {

	taus->Write(outFile);
	weights->Write(outFile);
	omegas->Write(outFile);
	sensorWeights->Write(outFile);

	fingerLengths->Write(outFile);
	fingerRadii->Write(outFile);
	fingerSpacings->Write(outFile);
	innerFingerRanges->Write(outFile);

	(*outFile) << coreDist		<< "\n";
	(*outFile) << leftShoulderDist	<< "\n";
	(*outFile) << rightShoulderDist	<< "\n";

	(*outFile) << fitness		<< "\n";
	
	(*outFile) << dominated		<< "\n";

	(*outFile) << validEval		<< "\n";
	(*outFile) << objectsEvaluated	<< "\n";
}

void GENOME::SendAsBrain(int fileIndex, int currentSlot) {

	char tempFileName[200];
	if ( fileIndex==-1 )
		sprintf(tempFileName,"/tmp/Files%d_%d/temp.dat",simParams->randSeed,currentSlot);
	else
		sprintf(tempFileName,"../Data/%d/temp%d.dat",simParams->runIndex,simParams->randSeed);

	char fileName[300];

	char tausFileName[200];	
	sprintf(tausFileName,"/tmp/Files%d_%d/%s",simParams->randSeed,currentSlot,TAUS_OUT_FILENAME);
	ofstream *outFile = new ofstream(tempFileName);
	(*outFile) << "1\n";
	taus->Write(outFile);
	outFile->close();
	delete outFile;
	if ( fileIndex==-1 )
		simParams->FileRename(	tempFileName,
					tausFileName);
	else {
		sprintf(fileName,"../Data/%d/%d_Taus_%d.dat",simParams->runIndex,simParams->randSeed,fileIndex);
		simParams->FileRename(tempFileName,fileName);
	}

	char weightsFileName[200];
	sprintf(weightsFileName,"/tmp/Files%d_%d/%s",simParams->randSeed,currentSlot,WEIGHTS_OUT_FILENAME);
	outFile = new ofstream(tempFileName);
	(*outFile) << "1\n";
	weights->Write(outFile);
	outFile->close();
	delete outFile;
	if ( fileIndex==-1 )
		simParams->FileRename(tempFileName,
					weightsFileName);
	else {
		sprintf(fileName,"../Data/%d/%d_Weights_%d.dat",simParams->runIndex,simParams->randSeed,fileIndex);
		simParams->FileRename(tempFileName,fileName);
	}

	char omegasFileName[200];
	sprintf(omegasFileName,"/tmp/Files%d_%d/%s",simParams->randSeed,currentSlot,OMEGAS_OUT_FILENAME);
	outFile = new ofstream(tempFileName);
	(*outFile) << "1\n";
	omegas->Write(outFile);
	outFile->close();
	delete outFile;
	if ( fileIndex==-1 )
		simParams->FileRename(tempFileName,
					omegasFileName);
	else {
		sprintf(fileName,"../Data/%d/%d_Omegas_%d.dat",simParams->runIndex,simParams->randSeed,fileIndex);
		simParams->FileRename(tempFileName,fileName);
	}


	char sensorsFileName[200];
	sprintf(sensorsFileName,"/tmp/Files%d_%d/%s",simParams->randSeed,currentSlot,SENSORS_OUT_FILENAME);

	outFile = new ofstream(tempFileName);
	(*outFile) << "1\n";
	sensorWeights->Write(outFile);
	outFile->close();
	delete outFile;
	if ( fileIndex==-1 )
		simParams->FileRename(tempFileName,
					sensorsFileName);
	else {
		sprintf(fileName,"../Data/%d/%d_Sensors_%d.dat",simParams->runIndex,simParams->randSeed,fileIndex);
		simParams->FileRename(tempFileName,fileName);
	}
}

void GENOME::SendEarlyStoppingData(int currentSlot) {

	char tempFileName[200];

	char fileName[300];
	ofstream *outFile;

	sprintf(tempFileName,"/tmp/Files%d_%d/temp1.dat",simParams->randSeed,currentSlot);
	outFile = new ofstream(tempFileName);
	objectOrdering->Write(outFile);
	outFile->close();
	delete outFile;
	outFile = NULL;
	sprintf(fileName,"/tmp/Files%d_%d/objectOrdering.dat",simParams->randSeed,currentSlot);
	simParams->FileRename(tempFileName,fileName);


	sprintf(tempFileName,"/tmp/Files%d_%d/temp3.dat",simParams->randSeed,currentSlot);
	outFile = new ofstream(tempFileName);
	coreDists->Write(outFile);
	outFile->close();
	delete outFile;
	outFile = NULL;
	sprintf(fileName,"/tmp/Files%d_%d/myCoreDists.dat",simParams->randSeed,currentSlot);
	simParams->FileRename(tempFileName,fileName);

	sprintf(tempFileName,"/tmp/Files%d_%d/temp4.dat",simParams->randSeed,currentSlot);
	outFile = new ofstream(tempFileName);
	leftShoulderDists->Write(outFile);
	outFile->close();
	delete outFile;
	outFile = NULL;
	sprintf(fileName,"/tmp/Files%d_%d/myLeftShoulderDists.dat",simParams->randSeed,currentSlot);
	simParams->FileRename(tempFileName,fileName);

	sprintf(tempFileName,"/tmp/Files%d_%d/temp5.dat",simParams->randSeed,currentSlot);
	outFile = new ofstream(tempFileName);
	rightShoulderDists->Write(outFile);
	outFile->close();
	delete outFile;
	outFile = NULL;
	sprintf(fileName,"/tmp/Files%d_%d/myRightShoulderDists.dat",simParams->randSeed,currentSlot);
	simParams->FileRename(tempFileName,fileName);

	sprintf(tempFileName,"/tmp/Files%d_%d/temp6.dat",simParams->randSeed,currentSlot);
	outFile = new ofstream(tempFileName);
	(*outFile) << objectsEvaluated;
	outFile->close();
	delete outFile;
	outFile = NULL;
	sprintf(fileName,"/tmp/Files%d_%d/objectsEvaluated.dat",simParams->randSeed,currentSlot);
	simParams->FileRename(tempFileName,fileName);
}

int GENOME::Successful(void) {

	int successCoreDist		= 	int(round(1000.0*coreDist)) >= 
						int(round(1000.0*CORE_DIST_MAX));

	int successLeftShoulderDist	= 	int(round(1000.0*leftShoulderDist))    >= 
						int(round(1000.0*LEFT_SHOULDER_DIST_MAX));

	int successRightShoulderDist	= 	int(round(1000.0*rightShoulderDist))  >= 
						int(round(1000.0*RIGHT_SHOULDER_DIST_MAX));

	return( successCoreDist && successLeftShoulderDist && successRightShoulderDist );
}

int  GENOME::SuccessfulOnAllEarlierObjects(void) {

	if ( objectsEvaluated == 1 )
		return( true );

	else {
		int failureFound = false;
		int currentObject = 0;

		while ( (!failureFound) && (currentObject<(objectsEvaluated-1)) ) {

			if ( SuccessfulOnObject(currentObject) )
				currentObject++;
			else
				failureFound = true;
		}

		return( !failureFound );
	}
}

int GENOME::SuccessfulOnObject(int objectIndex) {

	int successCoreDist		=	int(round(1000.0*coreDists->Get(0,objectIndex))) >= 
						int(round(1000.0*CORE_DIST_MAX));

	int successLeftShoulderDist	= 	int(round(1000.0*leftShoulderDists->Get(0,objectIndex)))    >= 
						int(round(1000.0*LEFT_SHOULDER_DIST_MAX));

	int successRightShoulderDist	= 	int(round(1000.0*rightShoulderDists->Get(0,objectIndex)))  >= 
						int(round(1000.0*RIGHT_SHOULDER_DIST_MAX));

	return( successCoreDist && successLeftShoulderDist && successRightShoulderDist );
}

void GENOME::SwapParentPointer(GENOME *firstParent, GENOME *secondParent) {

	if ( myParent ) {
		if ( myParent==firstParent )
			myParent = secondParent;
		else if ( myParent==secondParent )
			myParent = firstParent;
	}
}

int  GENOME::Superior(GENOME *other) {

	int lessExperience 	= objectsEvaluated < other->objectsEvaluated;

	if ( lessExperience )
		return( false );

	return( Dominates(other) );
}

int GENOME::SuperiorTo(GENOME *other) {

	Superior(other);
/*
        int lessExperience       = objectsEvaluated < other->objectsEvaluated;

        if ( lessExperience )
                return( false );

        int superiorForCore =   int(round(1000.0*       coreDists->Mean(0,0,0,objectsEvaluated-1))) >
                                int(round(1000.0*other->coreDists->Mean(0,0,0,other->objectsEvaluated-1)));

        int superiorForLeft =   int(round(1000.0*       leftShoulderDists->Mean(0,0,0,objectsEvaluated-1))) >
                                int(round(1000.0*other->leftShoulderDists->Mean(0,0,0,other->objectsEvaluated-1)));

        int superiorForRight =  int(round(1000.0*       rightShoulderDists->Mean(0,0,0,objectsEvaluated-1))) >
                                int(round(1000.0*other->rightShoulderDists->Mean(0,0,0,other->objectsEvaluated-1)));

        return( superiorForCore && superiorForLeft && superiorForRight );
*/

/*
	int minimumEvaluatedObjects = 1000;

	if ( objectsEvaluated < minimumEvaluatedObjects )

		minimumEvaluatedObjects = objectsEvaluated;

	if ( other->objectsEvaluated < minimumEvaluatedObjects )

		minimumEvaluatedObjects = objectsEvaluated;

	int superiorForCore =	int(round(1000.0*       coreDists->Mean(0,0,0,minimumEvaluatedObjects-1))) >
				int(round(1000.0*other->coreDists->Mean(0,0,0,minimumEvaluatedObjects-1)));

	int superiorForLeft =	int(round(1000.0*       leftShoulderDists->Mean(0,0,0,minimumEvaluatedObjects-1))) >
				int(round(1000.0*other->leftShoulderDists->Mean(0,0,0,minimumEvaluatedObjects-1)));

	int superiorForRight =	int(round(1000.0*       rightShoulderDists->Mean(0,0,0,minimumEvaluatedObjects-1))) >
				int(round(1000.0*other->rightShoulderDists->Mean(0,0,0,minimumEvaluatedObjects-1)));

	return( superiorForCore && superiorForLeft && superiorForRight );
*/
	return 0; //RISI
}

void GENOME::UpdateParentsChildDeathDistribution(void) {

	if ( myParent ) {		
		myParent->childDeathDistribution->Add(0,int(objectOrdering->Get(0,objectsEvaluated-1)),1.0);
		myParent = NULL;
	}
}

int GENOME::ValidParent(void) {

	return( !dominated );
}

void GENOME::Write(void) {

	printf("%d %3.3f %3.3f %3.3f ",dominated,coreDist,leftShoulderDist,rightShoulderDist);
}

// ----------------------------- Private methods ---------------------------------

void GENOME::CalculateCoreDist(int objectIndex, int currentSlot) {

	char fileName[200];
	ifstream *inFile;

	sprintf(fileName,"/tmp/Files%d_%d/myCoreDists.dat",simParams->randSeed,currentSlot);
	while ( !simParams->FileExists(fileName) );
	inFile = new ifstream(fileName);
	delete coreDists;
	coreDists = new MATRIX(inFile);
	inFile->close();
	delete inFile;
	inFile = NULL;
	simParams->FileDelete(fileName);

	if ( simParams->selectForCoreDists )
		coreDist = coreDists->Mean();
	else
		coreDist = CORE_DIST_MAX;
}

void GENOME::CalculateLeftShoulderDist(int objectIndex, int currentSlot) {

	char fileName[200];
	ifstream *inFile;

	sprintf(fileName,"/tmp/Files%d_%d/myLeftShoulderDists.dat",simParams->randSeed,currentSlot);
	while ( !simParams->FileExists(fileName) );
	inFile = new ifstream(fileName);
	delete leftShoulderDists;
	leftShoulderDists = new MATRIX(inFile);
	inFile->close();
	delete inFile;
	inFile = NULL;
	simParams->FileDelete(fileName);

	if ( simParams->selectForLeftShoulderDists )
		leftShoulderDist = leftShoulderDists->Mean();
	else
		leftShoulderDist = LEFT_SHOULDER_DIST_MAX;
}

void GENOME::CalculateRightShoulderDist(int objectIndex, int currentSlot) {

	char fileName[200];
	ifstream *inFile;

	sprintf(fileName,"/tmp/Files%d_%d/myRightShoulderDists.dat",simParams->randSeed,currentSlot);
	while ( !simParams->FileExists(fileName) );
	inFile = new ifstream(fileName);
	delete rightShoulderDists;
	rightShoulderDists = new MATRIX(inFile);
	inFile->close();
	delete inFile;
	inFile = NULL;
	simParams->FileDelete(fileName);

	if ( simParams->selectForRightShoulderDists )
		rightShoulderDist = rightShoulderDists->Mean();
	else
		rightShoulderDist = RIGHT_SHOULDER_DIST_MAX;
}


void GENOME::Clear(void) {

	if ( evalLengths ) {
		delete evalLengths;
		evalLengths = NULL;
	}

	if ( coreDists ) {
		delete coreDists;
		coreDists = NULL;
	}

	if ( leftShoulderDists ) {
		delete leftShoulderDists;
		leftShoulderDists = NULL;
	}

	if ( rightShoulderDists ) {
		delete rightShoulderDists;
		rightShoulderDists = NULL;
	}

	if ( objectOrdering ) {
		delete objectOrdering;
		objectOrdering = NULL;
	}
}

void GENOME::Initialize(void) {

	coreDists 		= new MATRIX(1,objectsToEvaluate,CORE_DIST_MAX);
	coreDist		= 0.0;

	leftShoulderDists 	= new MATRIX(1,objectsToEvaluate,LEFT_SHOULDER_DIST_MAX);
	leftShoulderDist	= 0.0;

	rightShoulderDists 	= new MATRIX(1,objectsToEvaluate,RIGHT_SHOULDER_DIST_MAX);
	rightShoulderDist	= 0.0;

	dominated		= false;
	validEval		= true;
	objectsEvaluated	= 0;

	objectOrdering 		= new MATRIX(1,objectsToEvaluate,0.0);
	for (int j=0;j<objectsToEvaluate;j++)
		//objectOrdering->Set(0,j,simParams->objectOrder->Get(0,(objectsToEvaluate-1)-j));
		objectOrdering->Set(0,j,simParams->objectOrder->Get(0,j));
		
	evalLengths = new MATRIX(1,objectsToEvaluate,0.0);
}

void GENOME::Invalidate(void) {

	coreDist = 0.0;
	leftShoulderDist = 0.0;
	rightShoulderDist = 0.0;
	coreDists->ReZero();
	leftShoulderDists->ReZero();
	rightShoulderDists->ReZero();

	validEval = false;
	fitness = coreDist * leftShoulderDist * rightShoulderDist;
}

#endif
