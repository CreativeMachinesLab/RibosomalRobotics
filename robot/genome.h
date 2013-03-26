/* ---------------------------------------------------
   FILE:     bodyPlan.h
	AUTHOR:   Josh Bongard
	DATE:     October 2, 2000
	FUNCTION: This class contains all information for
				 a single physical segment of an organism
				 in the MathEngine environment.
 -------------------------------------------------- */

#ifndef _GENOME_H
#define _GENOME_H

#include "matrix.h"

class GENOME {

public:

	MATRIX *taus;
	MATRIX *weights;
	MATRIX *omegas;
	MATRIX *sensorWeights;

	int	totalNeurons;

	MATRIX *fingerLengths;
	MATRIX *fingerRadii;
	MATRIX *fingerSpacings;
	MATRIX *innerFingerRanges;

	MATRIX *coreDists;
	double coreDist;

	MATRIX *leftShoulderDists;
	double leftShoulderDist;

	MATRIX *rightShoulderDists;
	double rightShoulderDist;

	double fitness;
	
	int dominated;

	int validEval;
	int objectsEvaluated;
	int objectsToEvaluate;

	MATRIX *objectOrdering;

	MATRIX *childDeathDistribution;

	GENOME *myParent;
	MATRIX *evalLengths;

	int    evalsCreated;
	int    evalsAtMoveUp;

public:

	GENOME(int evalsCurrent);
	GENOME(GENOME *parent);
	~GENOME(void);
	int  Age(int evalsCurrent);
	int  AllObjectsEvaluated(void);
	void Balance(int objectsOnLayer);
	void CalculateFitness(int currentGeneration, int rowIndex,
				int currentSlot);
	void Compress(void);
	int  Dominates(GENOME *otherGenome);
	void Expand(void);
	double MeanEvalTimes(void);
	void Mutate(void);
	void Print(void);
	void Reset(void);
	void Save(ofstream *outFile);
	void SendAsBrain(int fileIndex, int currentSlot);
	void SendEarlyStoppingData(int currentSlot);
	int  Successful(void);
	int  SuccessfulOnAllEarlierObjects(void);
	int  SuccessfulOnObject(int objectIndex);
	int  Superior(GENOME *other);
	int  SuperiorTo(GENOME *other);
	void SwapParentPointer(GENOME *firstParent, GENOME *secondParent);
	void UpdateParentsChildDeathDistribution(void);
	int  ValidParent(void);
	void Write(void);
	double round(double f) { return floor(f);}

private:
	void Clear(void);
	void CalculateCoreDist(int objectIndex, 	int currentSlot);
	void CalculateLeftShoulderDist(int objectIndex, int currentSlot);
	void CalculateRightShoulderDist(int objectIndex,int currentSlot);
	void Initialize(void);
	void Invalidate(void);
};

#endif
