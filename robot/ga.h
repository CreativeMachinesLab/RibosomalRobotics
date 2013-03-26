/* ---------------------------------------------------
   FILE:     growGA.h
	AUTHOR:   Josh Bongard
	DATE:     October 20, 2000
	FUNCTION: This class contains all information for
			  a population of variable-length genotypes

 -------------------------------------------------- */

#ifndef _GA_H
#define _GA_H

#include "iostream"
#include "fstream"
#include "genome.h"

class GA {

public:
	int      currGeneration;
	int	 currGenome;
	int      totalGenerations;
	GENOME   **genomes;
	char	 dataFileName[200];
	int	 numReplacements;
	int	 gensSinceLastReplacement;
	double	 meanObjectsEvaluated;
	double	 meanEvalTimes;
	double	 meanCoreDists;
	double	 meanLeftShoulderDists;
	double   meanRightShoulderDists;
	double	 bestCoreDist;
	double	 bestLeftShoulderDist;
	double	 bestRightShoulderDist;
	double   numGenomesEvaluated;
	MATRIX	 *dataForOutputFile;
	int	 rowIndex;
	int	 newNonDominatedFound;
	double   currBestFit;
	double	 previousBestFit;
	int      outputDataWritten;
	double   nums;
	int	 reinitializationMode;
	int	 reinitializationIndividual;
	int	 *ageLayerOfGenomeInTempSlots;

private:
	int	*slotStates;
	int	*genomeIndices;
	int	*objectsForLayer;

	MATRIX 	*stance_Prone;
	MATRIX 	*stance_Upright;
	MATRIX  *stance_Squat;
	MATRIX 	*stance_Squat_Upright;
	MATRIX  *stance_Four_Feet;
	MATRIX  *stance;
	MATRIX  *endStance;
	double	endStanceTransitionWhen;
	double  startRetraction;
	double  endExtension;
        double ratio;
        double endRatio;

public:

	GA(int loadState);
	~GA(void);

	void	Evolve();
	void	SaveState(void);
	void	SendBody(	int fileIndex, 
				int objectIndex,
				int currentSlot,
				int objsOnMyLayer);

private:
	int	AllSlotsEmpty(void);
	int     AgeLayer(int genomeIndex);
	int     AllDominated(void);
	int	AllEmpty(void);
	int     AllFinished(int firstSlot, int lastSlot);
	int     AllTooOld(int layerIndex);
	void    AttemptReplacement(int victimLayer, int aggressor, int slotIndex);
//	void	BalanceLayer(int genomeIndex, int layerIndex, int slotIndex);
	void    BalanceLayer(int layerIndex, int slotIndex);
	void    BalanceLayers(int layerIndex, int slotIndex);
	int     BestGenomeOnLayer(int ageLayer);
	int     BottomLayer(int layerIndex);
	void    CreateAndEvaluateChild(int parentIndex, int ageLayerOfChild, int slotIndex);
	void    CreateNewRandomGenome(int slotIndex);
	void    Discard(int genomeIndex);
	int     DiscardUnbalancedIncomingGenome(int aggressor, int victimLayer);
	int     Dominated(GENOME *currGenome, int layerIndex);
	int     Empty(int genomeIndex);
	void	EmptyAllSlots(void);
	int     EmptySlot(int currentSlot);
	void    EmptySlotOnce(int currentSlot);
	void    EvaluateAllLayers(void);
	void    EvaluateGenome(int genomeIndex, int slotIndex);
	void    EvaluateGenomeOnce(int genomeIndex, int slotIndex);
	void    EvaluateLayer(int layerIndex, int skipDominated, char *evaluatedNote);
	void	ExpandForAscension(int i, int slotIndex);
	void    ExpandGenome(int i, int slotIndex);
	void    ExpandGenomesOnLayer(int layerIndex);
//	void    ExpandLayer(int layerIndex, int slotIndex);
	void    ExpandLayer(int layerIndex, int slotIndex);
	void	ExpandLayers(int layerIndex, int slotIndex);
	void    ExpandThisGenome(int genomeIndex, int layerIndex, int slotIndex);
	void    ExpandThisLayer(int genomeIndex, int layerIndex, int slotIndex);
	void	FillSlot(int currentSlot);
	void    FillSlot(int genomeToSend, int currentSlot);
	void    FillSlot(int genomeToSend, int currentSlot, int whichObject);
	int	FindAnyone(int ageLayer);
	int	FindBest(void);
	int     FindBestOnLayer(int layerIndex);
	int	FindDominated(int ageLayer);
	int     FindGenomeOnActivatedLayer(void);
	int     FindHole(int i);
//	int	FindInferior(int i, int ageLayer);
	int     FindInferior(int aggressor, int layerIndex);
	int  	FindLayerForNewChild(void);	
	int	FindParent(int layerIndex);
	int	FindNonDominatedAndNotTooOld(int layerIndex);
	int	FindNonDominatedAndTooOld(int layerIndex);
	int     FindSlot(void);
	int	FindSuccessful(int layerIndex);
	int     FindSuccessfulGenome(void);
	int	FindTarget(int layerIndex);
	int	FindTooOld(int layerIndex);
	int	FindTooOldGenome(void);
	int     FindVictim(int i, int layerIndex);
	void 	FixLayer(int layerIndex, int slotIndex);
	void    FixLayers(int slotIndex);
	int	GenomesToEvaluateOnThisLayer(int layerIndex);
	double  GetFingerAngle(int fingerIndex, MATRIX *fingerSpacings);
	int     GetSlotState(int currentSlot);
	int	HandleAggressorDominatingVictim(int aggressor, int victimIndex, int slotIndex);
	int	HandleNeitherDominatingOneAnother(int aggressor, int victimIndex);
	void    HandleSuccess(void);
	void	HandleSuccessfulGenome(int successfulIndex, int slotIndex);
	int	InferiorOnLayer(int genomeIndex, int layerIndex);
	int	InTopLayer(int genomeIndex);
	int     InvalidGenomeOnLayer(int layerIndex);
	int     IsAVictim(int potentialVictim, int aggressor);
	int     IsNonDominatedInLayer(int genomeIndex, int ageLayer);
	void    KillDominated(void);
	int     LayerActivated(int layerIndex);
	int	LayerEvaluated(int layerIndex);
	int     MaxAgeForLayer(int layerIndex);
	int	MaxEvalsOnLayer(int layerIndex);
	double  MaxMult(void);
	int     MaxMultIndex(void);
	int	MaxSuccessfulObjectsOnLayer(int layerIndex);
	int	MostObjectsBeingTackled(void);
	void    Move(int from, int to);
	int	MovedLongEnoughAgo(int genomeIndex);
	void    MoveReinitializationPointer(void);
	int	NoSuccessesForObject(int objectIndex, int layerIndex);
	int     NotEnoughExperienceForLayer(int genomeIndex, int layerIndex);
	int     NumFilled(int ageLayer);
	int     NumNonDominated(int ageLayer);
	int     NumNonDominated(GENOME *currGenome, int ageLayer);
	int	ObjectsForLayer(int layerIndex);
	void    PrintGenome(int genomeIndex, char *note);
	void	PrintLayer(int layerIndex);
	void	PrintSlots(void);
	int     RandomGenome(int ageLayer);
	void    RandomizeGenome(int genomeIndex);
	void    RecalculateFront(int ageLayer);
	void    RecalculateFrontStrict(int ageLayer);
	void    RecalculateFrontsStrict(int victimLayer);
	void	RecalculateFrontUsingDomination(int ageLayer);
	int     ReinitializationFinished(void);
	void    ReplaceVictimAbove(int i, int slotIndex);
	void	ResetAllGenomes(void);
	void    SaveBest(void);
	void    SaveObjsConquered(void);
	void    SaveOutData(void);
	void    SendArm_Joints(ofstream *outFile);
	void    SendArms(ofstream *outFile);
	void	SendCore(ofstream *outFile);
	void	SendEarlyStoppingData(GENOME *currGenome, int currentSlot);
	void	SendHipGirdles(ofstream *outFile);
	void    SendHips(ofstream *outFile);
	void    SendLeg_Joints(ofstream *outFile);
	void    SendLegs(ofstream *outFile);
	void	SendLowerArm(ofstream *outFile);
	void	SendLowerArm_Joint(ofstream *outFile);
	void    SendPostscript(ofstream *outFile);
	void    SendSegments(ofstream *outFile);
	void    SendShoulder(ofstream *outFile);
	void	SendShoulderGirdles(ofstream *outFile);
	void    SendShoulders(ofstream *outFile);
	void    SendSpine_Joints(ofstream *outFile);
	void	SendTargetObject(ofstream *outFile, int objectIndex);
	void	SendUpperArm(ofstream *outFile);
	void	SendUpperArm_Joint(ofstream *outFile);
	void    SeverChildConnections(GENOME *dyingParent);
	int     SlotFinished(int currentSlot);
	void    SortGenomesOnLayer(int layerIndex);
	void	SortGenomesOnLayers(void);
	void	Stance_Create(void);
	void	Stance_Destroy(void);
	void    Stance_SetForBetweenScaffolding(int objectIndex);
	void    Stance_SetForNoScaffolding(void);
	void    Stance_SetForWithinScaffolding(int objectIndex);
	void    Stance_SetForWithinScaffolding2(int objectIndex);
	void    Stance_SetForWithinScaffolding3(int objectIndex);
	void    Stance_SetForWithinScaffolding4(int objectIndex);
	int	SuccessfulGenomePresent(void);
	int	SuperiorOnLayer(int genomeIndex, int layerIndex);
	void    SwapParentPointers(GENOME *firstParent, GENOME *secondParent);
	int	TooMuchExperienceForLayer(int genomeIndex, int layerIndex);
	int     TooOld(int genomeIndex);
	void    TryMoveUp(int i, int slotIndex);
	void    TurnOffReinitialization(void);
	void    TurnOnReinitialization(void);
	int     VictimFailed(int aggressor, int victimIndex, int slotIndex);
	void	WaitForSensorReadings(void);
	void    WriteBest(ofstream *outFile, int layerIndex);
	void    WriteMean(ofstream *outFile, int layerIndex);
	void    WriteAll(void);
	int     YoungEnoughParent(int genomeIndex);
};

#endif
