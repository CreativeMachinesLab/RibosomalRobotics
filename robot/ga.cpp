/* ---------------------------------------------------
   FILE:     growGA.cpp
	AUTHOR:   Josh Bongard
	DATE:     October 20, 2000
	FUNCTION: This class contains all information for
				 a population of variable-length genotypes
 -------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#ifndef _GA_CPP
#define _GA_CPP

#include "constants.h"
#include "genome.h"
#include "ga.h"
#include "simParams.h"

extern SIM_PARAMS	*simParams;
extern char		TEMP_FILENAME[100];

GA::GA(int loadState) {

	Stance_Create();

	ifstream *inFile = NULL;
	if ( loadState ) {
		char fileName[200];
		inFile = new ifstream(fileName);
	}

	slotStates    = new int[simParams->numSlots];
	genomeIndices = new int[simParams->numSlots];

	char dirName[200];
	for (int i=0;i<simParams->numSlots;i++) {
		sprintf(dirName,"%s%d_%d",TEMP_FILENAME_PREFIX_NEW,simParams->randSeed,i);
		simParams->DirectoryMake(dirName);
		slotStates[i]    = EMPTY;
		genomeIndices[i] = -1;
	}	

	sprintf(dataFileName,"../Data/%d/runData_%d.dat",simParams->runIndex,simParams->randSeed);
	if ( !loadState ) {
		simParams->FileDelete(dataFileName);
		ofstream *outFile = new ofstream(dataFileName);
		outFile->close();
		delete outFile;
	}

/*
	sprintf(simParams->bodyParamsFileName,"../Data/%d/bodyParams_%d.dat",simParams->runIndex,simParams->randSeed);
	if ( !loadState ) {
		simParams->FileDelete(simParams->bodyParamsFileName);
		ofstream *outFile = new ofstream(simParams->bodyParamsFileName);
		outFile->close();
		delete outFile;
	}

	sprintf(simParams->objsConqFileName,"../Data/%d/objsConq_%d.dat",simParams->runIndex,simParams->randSeed);
	if ( !loadState ) {
		simParams->FileDelete(simParams->objsConqFileName);
		ofstream *outFile = new ofstream(simParams->objsConqFileName);
		outFile->close();
		delete outFile;
	}
*/
	genomes = new GENOME * [simParams->popSize + simParams->numSlots];
	for (int i=0;i<(simParams->popSize+simParams->numSlots);i++)
		if ( i<simParams->popSize )
			genomes[i] = new GENOME(int(numGenomesEvaluated));
		else
			genomes[i] = NULL;

	dataForOutputFile = new MATRIX(DATA_FILE_BUFFER,simParams->dataFileBufferWidth,0.0);
	rowIndex = 0;

	numGenomesEvaluated = 0.0;
	outputDataWritten   = false;

	meanObjectsEvaluated 	= 0.0;
	meanEvalTimes	     	= 0.0;
	meanCoreDists	     	= 0.0;
	meanLeftShoulderDists			= 0.0;
	meanRightShoulderDists		= 0.0;
	bestCoreDist		= 0.0;
	bestLeftShoulderDist			= 0.0;
	bestRightShoulderDist		= 0.0;
	nums		     	= 0.0;
	ageLayerOfGenomeInTempSlots	= new int[simParams->numSlots];
	previousBestFit = 0.0;
	currBestFit     = 0.0;

	reinitializationMode = false;

	if ( loadState ) {

		(*inFile) >> currGeneration;

		inFile->close();
		delete inFile;
		inFile = NULL;
	}

	objectsForLayer = new int[simParams->numLayers];
	for (int i=0;i<simParams->numLayers;i++)
		objectsForLayer[i] = 1;
}

GA::~GA(void) {

	Stance_Destroy();

	delete[] ageLayerOfGenomeInTempSlots;
	ageLayerOfGenomeInTempSlots = NULL;

	if ( objectsForLayer ) {
		delete objectsForLayer;
		objectsForLayer = NULL;
	}

	if ( slotStates ) {
		delete slotStates;
		slotStates = NULL;
	}

	if ( genomeIndices ) {
		delete genomeIndices;
		genomeIndices = NULL;
	}

	if ( dataForOutputFile ) {
		delete dataForOutputFile;
		dataForOutputFile = NULL;
	}

	for (int i=0;i<(simParams->popSize+simParams->numSlots);i++)
		if ( genomes[i] )
			genomes[i]->myParent = NULL;

	for (int i=0;i<(simParams->popSize+simParams->numSlots);i++) {
		if ( genomes[i] ) {
			delete genomes[i];
			genomes[i] = NULL;
		}
	}
	delete[] genomes;
	genomes = NULL;
}	

void GA::Evolve(void) {
	
	EvaluateAllLayers();

	int currentSlot = 0;
	while ( 1 ) {

		if ( GetSlotState(currentSlot)==EMPTY ) {
                	int i = FindLayerForNewChild();
			if ( BottomLayer(i) && reinitializationMode && (simParams->numLayers>1)) {
				CreateNewRandomGenome(currentSlot);
				printf("\n");
			}
                       	else if ( BottomLayer(i) && AllTooOld(i) && (simParams->numLayers>1) ) {
                               	TurnOnReinitialization();
                               	CreateNewRandomGenome(currentSlot);
				printf("\n");
                       	}
			else {
	                	int parentIndex = FindParent(i);
				if ( parentIndex > -1 )
	       		         	CreateAndEvaluateChild(parentIndex,i,currentSlot);
			}
		}
		else if ( GetSlotState(currentSlot)==FINISHED ) {
                        int aggressor = EmptySlot(currentSlot);
                        if ( GetSlotState(currentSlot)==EMPTY ) {
				AttemptReplacement(AgeLayer(aggressor),aggressor,currentSlot);
				printf("\n");
			}
		}
		currentSlot++;
		if ( currentSlot == simParams->numSlots )
			currentSlot=0;
	}
}

void GA::SaveState(void) {

	char fileName[200];
	sprintf(fileName,"../Data/%d/%d_savedState.dat",simParams->runIndex,simParams->randSeed);
	ofstream *outFile = new ofstream(fileName);

	(*outFile) << -1 	<< "\n";

	for (int g=0;g<simParams->popSize;g++)
		genomes[g]->Save(outFile);

	(*outFile) << currGeneration  		<< "\n";

	outFile->close();
	delete outFile;
	outFile = NULL;
}

void GA::SendBody(	int fileIndex, 
			int objectIndex,
			int currentSlot,
                        int objsOnMyLayer) {

        if ( simParams->betweenScaffolding )
                Stance_SetForBetweenScaffolding(simParams->stanceIndex);

        else if ( simParams->withinScaffolding )
                Stance_SetForWithinScaffolding(simParams->stanceIndex);

        else if ( simParams->withinScaffolding2 )
                Stance_SetForWithinScaffolding2(simParams->stanceIndex);

        else if ( simParams->withinScaffolding3 )
                Stance_SetForWithinScaffolding3(simParams->stanceIndex);

        else if ( simParams->withinScaffolding4 )
                Stance_SetForWithinScaffolding4(simParams->stanceIndex);
        else
                Stance_SetForNoScaffolding();

	char tempFileName[200];
	if ( fileIndex==-1 )
		sprintf(tempFileName,"/tmp/Files%d_%d/temp.dat",simParams->randSeed,currentSlot);
	else
		sprintf(tempFileName,"../Data/%d/temp%d.dat",simParams->runIndex,simParams->randSeed);

	ofstream *outFile = new ofstream(tempFileName);

	(*outFile) << "1\n";
	(*outFile) << 20 << " " << 18 << "\n";

	SendTargetObject(outFile,0);

	SendCore(outFile);

	SendSegments(outFile);

	SendShoulderGirdles(outFile);
	SendShoulders(outFile);
	SendArms(outFile);

	SendHipGirdles(outFile);
	SendHips(outFile);
	SendLegs(outFile);
	
	SendSpine_Joints(outFile);
	SendArm_Joints(outFile);
	SendLeg_Joints(outFile);

	SendPostscript(outFile);

	outFile->close();
	delete outFile;
	outFile = NULL;

	char bodyFileName[200];
	sprintf(bodyFileName,"/tmp/Files%d_%d/%s",simParams->randSeed,currentSlot,BODY_OUT_FILENAME);

	if ( fileIndex==-1 )
		simParams->FileRename(tempFileName,bodyFileName);
	else {
		char fileName[200];
		sprintf(fileName,"../Data/%d/%d_Body_%d.dat",simParams->runIndex,simParams->randSeed,fileIndex);
		simParams->FileRename(tempFileName,fileName);
	}
}

// ----------------------------------------------------------------
//                           Private methods
// ----------------------------------------------------------------

int  GA::AllSlotsEmpty(void) {

	int allEmpty = true;
	int currentSlot=0;

	while ( allEmpty && (currentSlot<simParams->numSlots) ) {

		if ( GetSlotState(currentSlot)==EMPTY )
			currentSlot++;
		else
			allEmpty = false;
	}

	return( allEmpty );
}

int  GA::AgeLayer(int genomeIndex) {

	if ( genomeIndex>=simParams->popSize )

		return( ageLayerOfGenomeInTempSlots[genomeIndex-simParams->popSize] );

	else {

		double rawFraction = double(genomeIndex)/double(simParams->numIndsPerLayer);
		double layer = int( floor(rawFraction) );

		if ( layer >= simParams->numLayers )
			layer = simParams->numLayers-1;

		return( int(layer) );
	}
}

int  GA::AllDominated(void) {

	int allDominated = true;

	int i=0;
	while ( (i<simParams->popSize) && (allDominated) ) {

		if ( genomes[i] )
			if ( !genomes[i]->dominated )
				allDominated = false;
		i++;
	}

	return( allDominated );
}

int  GA::AllEmpty(void) {

	int allEmpty = true;
	int i=0;

	while ( (i<simParams->popSize) && (allEmpty) ) {

		if ( GetSlotState(i)!=EMPTY )
			allEmpty = false;
		else
			i++;
	}

	return( allEmpty );
}

int GA::AllFinished(int firstSlot, int lastSlot) {

	int allFinished = true;

	int i=firstSlot;

	while ( allFinished && (i<=lastSlot) )
		if ( !SlotFinished(i) )
			allFinished = false;
		else
			i++;

	return( allFinished );
}

int  GA::AllTooOld(int layerIndex) {

	int allTooOld = true;

	for (int i=0;i<simParams->popSize;i++)

		if ( 	genomes[i] &&
			(AgeLayer(i)==layerIndex) &&
			(genomes[i]->Age(int(numGenomesEvaluated)) <= MaxAgeForLayer(layerIndex)) )

			allTooOld = false;
			
	return( allTooOld );
}

void GA::AttemptReplacement(int victimLayer, int aggressor, int slotIndex) {

	int discarded = DiscardUnbalancedIncomingGenome(aggressor,victimLayer);
	if ( discarded )
		return;

        int victimIndex = FindVictim(aggressor,victimLayer);
	if ( victimIndex==-1 ) {
		Discard(aggressor);
		return;
	}

        TryMoveUp(victimIndex,slotIndex);

	if ( !Empty(victimIndex) )

		Discard(victimIndex);

	Move(aggressor,victimIndex);

	HandleSuccessfulGenome(simParams->successfulGenomeIndex,slotIndex);

	BalanceLayers(victimLayer,slotIndex);

	ExpandLayers(victimLayer,slotIndex);

	RecalculateFrontsStrict(victimLayer);
}

void GA::BalanceLayer(int layerIndex, int slotIndex) {

	objectsForLayer[layerIndex] = MaxEvalsOnLayer(layerIndex);

	while (	(objectsForLayer[layerIndex]>1) &&
		(NoSuccessesForObject(objectsForLayer[layerIndex]-1,layerIndex)) &&
		(NoSuccessesForObject(objectsForLayer[layerIndex]-2,layerIndex)) )

		objectsForLayer[layerIndex] = objectsForLayer[layerIndex] - 1;

	for (int i=0;i<simParams->popSize;i++)

		if (	(genomes[i])	&&
			(AgeLayer(i)==layerIndex) )

			genomes[i]->Balance(ObjectsForLayer(layerIndex));

	RecalculateFrontStrict(layerIndex);			
}

void GA::BalanceLayers(int layerIndex, int slotIndex) {

	for (int i=layerIndex;i<simParams->numLayers;i++)

		BalanceLayer(i,slotIndex);
}

int  GA::BestGenomeOnLayer(int ageLayer) {

	double bestFitness = -1000.0;
	int    bestIndex = -1;

	for (int i=0;i<simParams->popSize;i++)

		if ( genomes[i] && (AgeLayer(i)==ageLayer) ) {

			int considerGenome = genomes[i]->fitness > bestFitness;

			if (	considerGenome ) { 
				bestFitness = genomes[i]->fitness;
				bestIndex = i;
			}
		}

	return( bestIndex );
}

int  GA::BottomLayer(int layerIndex) {

	return ( layerIndex==0 );
}

void GA::CreateAndEvaluateChild(int parentIndex, int ageLayerOfChild, int slotIndex) {

	int genomeIndex = simParams->popSize + slotIndex;
	genomes[genomeIndex] = new GENOME(genomes[parentIndex]);

	while ( genomes[genomeIndex]->objectsToEvaluate < ObjectsForLayer(ageLayerOfChild) )
		genomes[genomeIndex]->Expand();

	while ( genomes[genomeIndex]->objectsToEvaluate > ObjectsForLayer(ageLayerOfChild) )
		genomes[genomeIndex]->Compress();

	ageLayerOfGenomeInTempSlots[slotIndex] = ageLayerOfChild;

	genomes[genomeIndex]->Mutate();

	FillSlot(genomeIndex,slotIndex);
}

void GA::CreateNewRandomGenome(int slotIndex) {

	if ( !Empty(reinitializationIndividual) )

		TryMoveUp(reinitializationIndividual,slotIndex);

	if ( !Empty(reinitializationIndividual) )

		Discard(reinitializationIndividual);

	genomes[reinitializationIndividual] =

		new GENOME(int(numGenomesEvaluated));

        while ( genomes[reinitializationIndividual]->objectsToEvaluate < ObjectsForLayer(0) )

                genomes[reinitializationIndividual]->Expand();

	genomes[reinitializationIndividual]->Mutate();

	EvaluateGenome(reinitializationIndividual,slotIndex);

	char note[100];
	sprintf(note,"New random genome created.\n");
	PrintGenome(reinitializationIndividual,note);

	HandleSuccessfulGenome(reinitializationIndividual,slotIndex);

	BalanceLayers(0,slotIndex);

	ExpandLayers(0,slotIndex);

        RecalculateFrontsStrict(0);
/*
	FixLayer(AgeLayer(reinitializationIndividual),slotIndex);

	char note[100];
	sprintf(note,"New random genome created.");
	PrintGenome(reinitializationIndividual,note);
*/
	MoveReinitializationPointer();
	if ( ReinitializationFinished() )
		TurnOffReinitialization();
}

void GA::Discard(int genomeIndex) {

	char note[100];
	sprintf(note,"Discarding genome.");
	PrintGenome(genomeIndex,note);	

	delete genomes[genomeIndex];
	genomes[genomeIndex] = NULL;

	if ( simParams->successfulGenomeIndex == genomeIndex ) {
		simParams->successOccurred = false;
		simParams->successfulGenomeIndex = -1;
	}
}

int GA::DiscardUnbalancedIncomingGenome(int aggressor, int victimLayer) {

	if ( genomes[aggressor]->objectsToEvaluate != ObjectsForLayer(victimLayer) ) {

		Discard(aggressor);
		return ( true );
	}

	return( false );
}

int  GA::Dominated(GENOME *currGenome, int ageLayer) {

	int i=0;

	int dominated = false;

	while ( (i<simParams->popSize) && (!dominated) ) {

		if ( 	(genomes[i]) && 
			(AgeLayer(i)==ageLayer) &&
			(genomes[i]!=currGenome) ) {

			if ( genomes[i]->Dominates(currGenome) )
				dominated = true;
			else
				i++;
		}
		else
			i++;
	}

	return( dominated );
}

int  GA::Empty(int genomeIndex) {

	return( genomes[genomeIndex] == NULL );
}

void GA::EmptyAllSlots(void) {

	int currentSlot=0;
	while ( !AllSlotsEmpty() ) {
		if ( GetSlotState(currentSlot)==FINISHED ) {
        		int tmp = EmptySlot(currentSlot);
			if ( GetSlotState(currentSlot)==EMPTY )
				Discard(tmp);
		}
		currentSlot++;
		if ( currentSlot==simParams->numSlots )
			currentSlot=0;
	}
}

int GA::EmptySlot(int currentSlot) {

	int genomeIndex = genomeIndices[currentSlot];

	GENOME *currentGenome = genomes[genomeIndex];

	currentGenome->CalculateFitness(0,rowIndex,currentSlot);
	currentGenome->dominated = Dominated(currentGenome,AgeLayer(genomeIndex));

	int validEval = currentGenome->validEval;
	int moreObjs  = !currentGenome->AllObjectsEvaluated();
	int nonDominated = !currentGenome->dominated;
	int reload = validEval && moreObjs;

	if ( simParams->useEarlyStopping )
		reload = reload && nonDominated; 

	if ( reload ) {

		SendEarlyStoppingData(currentGenome,currentSlot);
		currentGenome->SendEarlyStoppingData(currentSlot);
		currentGenome->SendAsBrain(-1,currentSlot);
		SendBody(-1,	int(currentGenome->objectOrdering->Get(0,currentGenome->objectsEvaluated)),
				currentSlot,currentGenome->objectsToEvaluate);
		slotStates[currentSlot] = RUNNING;
	}
	else {
		genomeIndices[currentSlot] 	= -1;
		slotStates[currentSlot] 	= EMPTY;

		if ( !currentGenome->validEval ) {
			currentGenome->coreDist			= 0.0;
			currentGenome->leftShoulderDist		= 0.0;
			currentGenome->rightShoulderDist	= 0.0;
			currentGenome->dominated		= true;
		}

		numGenomesEvaluated  = numGenomesEvaluated + 1;

		simParams->totalTimeStepsUsed = simParams->totalTimeStepsUsed + 
						long((currentGenome->MeanEvalTimes()*currentGenome->objectsEvaluated));

		if ( (currentGenome->MeanEvalTimes()>1.0) && (currentGenome->coreDist>0.0) ) {
			meanObjectsEvaluated 	= meanObjectsEvaluated 		+ currentGenome->objectsEvaluated;
			meanEvalTimes 	     	= meanEvalTimes			+ currentGenome->MeanEvalTimes();
			meanCoreDists		= meanCoreDists			+ currentGenome->coreDist;
			meanLeftShoulderDists	= meanLeftShoulderDists		+ currentGenome->leftShoulderDist;
			meanRightShoulderDists	= meanRightShoulderDists	+ currentGenome->rightShoulderDist;
			nums		     	= nums + 1;

			if ( (currentGenome->fitness) > (bestCoreDist * bestLeftShoulderDist * bestRightShoulderDist) ) {
				bestCoreDist		= currentGenome->coreDist;
				bestLeftShoulderDist	= currentGenome->leftShoulderDist;
				bestRightShoulderDist	= currentGenome->rightShoulderDist;
			}
		}

		if ( currentGenome->Successful() ) {
			simParams->successOccurred = true;
			simParams->successfulGenomeIndex = genomeIndex;
		}
	}

	currentGenome = NULL;

	if ( (int(numGenomesEvaluated)%simParams->popSize)==0 )
		WriteAll();

	return genomeIndex;
}

void GA::EmptySlotOnce(int currentSlot) {

	int genomeIndex = genomeIndices[currentSlot];

	GENOME *currentGenome = genomes[genomeIndex];

	currentGenome->CalculateFitness(0,rowIndex,currentSlot);
	currentGenome->dominated = Dominated(currentGenome,AgeLayer(genomeIndex));

	genomeIndices[currentSlot] 	= -1;
	slotStates[currentSlot] 	= EMPTY;

	if ( !currentGenome->validEval ) {
		currentGenome->coreDist			= 0.0;
		currentGenome->leftShoulderDist		= 0.0;
		currentGenome->rightShoulderDist	= 0.0;
		currentGenome->dominated		= true;
	}

	numGenomesEvaluated  = numGenomesEvaluated + 1;

	simParams->totalTimeStepsUsed = simParams->totalTimeStepsUsed + 
					long(currentGenome->MeanEvalTimes()*currentGenome->objectsEvaluated);

	if ( (currentGenome->MeanEvalTimes()>1.0) && (currentGenome->coreDist>0.0) ) {
		meanObjectsEvaluated 	= meanObjectsEvaluated 		+ currentGenome->objectsEvaluated;
		meanEvalTimes 	     	= meanEvalTimes			+ currentGenome->MeanEvalTimes();
		meanCoreDists		= meanCoreDists			+ currentGenome->coreDist;
		meanLeftShoulderDists	= meanLeftShoulderDists		+ currentGenome->leftShoulderDist;
		meanRightShoulderDists	= meanRightShoulderDists	+ currentGenome->rightShoulderDist;
		nums		     	= nums + 1;

		if ( (currentGenome->fitness) > (bestCoreDist * bestLeftShoulderDist * bestRightShoulderDist) ) {
			bestCoreDist		= currentGenome->coreDist;
			bestLeftShoulderDist	= currentGenome->leftShoulderDist;
			bestRightShoulderDist	= currentGenome->rightShoulderDist;
		}
	}

	currentGenome = NULL;

        if ( (int(numGenomesEvaluated)%simParams->popSize)==0 )
                WriteAll();
}


void GA::EvaluateAllLayers(void) {

	int layerToEvaluate = true;

	while ( layerToEvaluate ) {

		layerToEvaluate = false;

		int bestLayerIndex = -1;
		double bestFitness = -1000.0;
		int    bestGenomeOnLayer;
		double bestFitnessOnLayer = -1;

		for (int layerIndex=0;layerIndex<simParams->numLayers;layerIndex++) {

			if ( !LayerEvaluated(layerIndex) ) {

				layerToEvaluate = true;

				bestGenomeOnLayer = BestGenomeOnLayer(layerIndex);

				bestFitnessOnLayer = genomes[bestGenomeOnLayer]->fitness;

				if ( bestFitnessOnLayer > bestFitness ) {

					bestFitness = bestFitnessOnLayer;
					bestLayerIndex = layerIndex;
				}
			}
		}

		if ( layerToEvaluate ) {
			char note[100];
			sprintf(note,"Initial random genome.");
			EvaluateLayer(bestLayerIndex,false,note);
			HandleSuccessfulGenome(simParams->successfulGenomeIndex,0);
		}
	}
}

void GA::EvaluateGenome(int genomeIndex, int slotIndex) {

	int whichObject = int(genomes[genomeIndex]->objectOrdering->Get(0,genomes[genomeIndex]->objectsEvaluated));

	FillSlot(genomeIndex,slotIndex,whichObject);

	int slotState = GetSlotState(slotIndex);

	while ( slotState!=EMPTY ) {

		if ( slotState==FINISHED )
			EmptySlot(slotIndex);
	
		slotState = GetSlotState(slotIndex);
	}

        if ( genomes[genomeIndex]->Successful() ) {
                simParams->successOccurred = true;
		simParams->successfulGenomeIndex = genomeIndex;
	}
}

void GA::EvaluateGenomeOnce(int genomeIndex, int slotIndex) {

	FillSlot(genomeIndex,slotIndex);

	int slotState = GetSlotState(slotIndex);

	while ( slotState!=EMPTY ) {

		if ( slotState==FINISHED )
			EmptySlotOnce(slotIndex);
	
		slotState = GetSlotState(slotIndex);
	}

        if (	genomes[genomeIndex]->Successful() &&
		genomes[genomeIndex]->objectsEvaluated==ObjectsForLayer(AgeLayer(genomeIndex)) ) {
                simParams->successOccurred = true;
		simParams->successfulGenomeIndex = genomeIndex;
	}

        if ( 	(genomes[genomeIndex]->Successful()) && 
		(genomes[genomeIndex]->objectsEvaluated==simParams->totalObjects) ) {
		WriteAll();
                exit(0);
        }

}

void GA::EvaluateLayer(int layerIndex, int skipDominated, char *evaluatedNote) {

	int currentSlot = 0;

	int currGenome = 0;
	while ( AgeLayer(currGenome)!=layerIndex )
		currGenome++;

	while ( !( 	(AgeLayer(currGenome)!=layerIndex) || 
			(currGenome>=simParams->popSize) ) ) {

		if ( GetSlotState(currentSlot)==EMPTY ) {
			if ( 	(AgeLayer(currGenome)==layerIndex) && 
				(currGenome<simParams->popSize) ) {

				if ( 	( (skipDominated) && 
					  (Dominated(genomes[currGenome],layerIndex)) ) ||
					(genomes[currGenome]->objectsEvaluated==genomes[currGenome]->objectsToEvaluate) ) {
                                        char note[100];
                                        sprintf(note,"No change.");
                                        PrintGenome(currGenome,note);
				}
				else {
					FillSlot(currGenome,currentSlot);
				}
				currGenome++;
			}
		}
		else if ( GetSlotState(currentSlot)==FINISHED ) {
			int i = EmptySlot(currentSlot);
			if ( GetSlotState(currentSlot)==EMPTY ) {
                        	RecalculateFrontStrict(AgeLayer(i));
				PrintGenome(i,evaluatedNote);
			}
		}
		currentSlot = currentSlot + 1;
		if ( currentSlot == simParams->numSlots )
			currentSlot = 0;
	}

	while ( !AllSlotsEmpty() ) {
                if ( GetSlotState(currentSlot)==FINISHED ) {
                        int i = EmptySlot(currentSlot);
                        if ( GetSlotState(currentSlot)==EMPTY ) {
                                RecalculateFrontStrict(AgeLayer(i));
                                PrintGenome(i,evaluatedNote);
                        }
                }
                currentSlot = currentSlot + 1;
                if ( currentSlot == simParams->numSlots )
                        currentSlot = 0;
	}
}

void GA::ExpandForAscension(int i, int slotIndex) {

	int nonDominated 	= !Dominated(genomes[i],AgeLayer(i)+1);
	int inexperienced 	= NotEnoughExperienceForLayer(i,AgeLayer(i)+1);

	while ( nonDominated && inexperienced ) {
		ExpandGenome(i,slotIndex);
		nonDominated = !Dominated(genomes[i],AgeLayer(i)+1);
		inexperienced = NotEnoughExperienceForLayer(i,AgeLayer(i)+1);
	}
}

void GA::ExpandGenome(int i, int slotIndex) {

	int genomeIndex = simParams->popSize + slotIndex;

	int temp = ageLayerOfGenomeInTempSlots[slotIndex];
	GENOME *tempGenome = genomes[genomeIndex];
	genomes[genomeIndex] = NULL;

	ageLayerOfGenomeInTempSlots[slotIndex] = AgeLayer(i);

	Move(i,genomeIndex);

	genomes[genomeIndex]->Expand();

	if ( simParams->successfulGenomeIndex == genomeIndex &&
		(!genomes[genomeIndex]->Successful()) ) {
		simParams->successOccurred = false;
		simParams->successfulGenomeIndex = -1;
	}

	char note[100];

	if ( !Dominated(genomes[genomeIndex],ageLayerOfGenomeInTempSlots[slotIndex]) ) {

		EvaluateGenomeOnce(genomeIndex,slotIndex);
		sprintf(note,"Expanded genome.");
	}
	else
		sprintf(note,"no change.");

	Move(genomeIndex,i);

	genomes[genomeIndex] = tempGenome;
	tempGenome = NULL;
	ageLayerOfGenomeInTempSlots[slotIndex] = temp;

	FixLayer(AgeLayer(i),slotIndex);

	PrintGenome(i,note);
}

void GA::ExpandGenomesOnLayer(int layerIndex) {

	for (int i=0;i<simParams->popSize;i++)

		if (	genomes[i]	&&
			(AgeLayer(i)==layerIndex) )

			while ( genomes[i]->objectsToEvaluate < objectsForLayer[layerIndex] )
				genomes[i]->Expand();
}

void GA::ExpandLayer(int layerIndex, int slotIndex) {

	int evaluationsForThisLayer = GenomesToEvaluateOnThisLayer(layerIndex);

	for (int i=0;i<simParams->popSize;i++) {

		int wasEvaled = false;

		while (	(genomes[i])	&&
			(AgeLayer(i)==layerIndex) &&
			(genomes[i]->objectsEvaluated<genomes[i]->objectsToEvaluate)	&&
			(!Dominated(genomes[i],layerIndex)) ) {

			EvaluateGenomeOnce(i,slotIndex);
			wasEvaled = true;
                	RecalculateFrontStrict(layerIndex);

			char note[100];
			sprintf(note,"Genome evaluated once.");
			PrintGenome(i,note);
		}
	}
}

void GA::ExpandLayers(int layerIndex, int slotIndex) {

	for (int i=layerIndex;i<simParams->numLayers;i++)

		ExpandLayer(i,slotIndex);
}

void  GA::ExpandThisGenome(int genomeIndex, int layerIndex, int slotIndex) {

	int relativeExperience =
		genomes[genomeIndex]->objectsToEvaluate -
		ObjectsForLayer(layerIndex);

	while (	relativeExperience<0 ) {
		ExpandGenome(genomeIndex, slotIndex);
		relativeExperience =
		        genomes[genomeIndex]->objectsToEvaluate -
		        ObjectsForLayer(layerIndex);
	}

	int nonDominated = 	!Dominated(genomes[genomeIndex],layerIndex);
	int objectsRemaining = 	genomes[genomeIndex]->objectsEvaluated <
				ObjectsForLayer(AgeLayer(genomeIndex));
	int evalsRemaining = 	nonDominated && objectsRemaining;

	while ( evalsRemaining ) {
		EvaluateGenomeOnce(genomeIndex,slotIndex);
		char note[100];
		sprintf(note,"Expanded genome.");
		PrintGenome(genomeIndex,note);
		nonDominated = 		!Dominated(genomes[genomeIndex],layerIndex);
		objectsRemaining = 	genomes[genomeIndex]->objectsEvaluated <
					ObjectsForLayer(AgeLayer(genomeIndex));
		evalsRemaining = 	nonDominated && objectsRemaining;
	}

	RecalculateFrontStrict(layerIndex);
}

void  GA::ExpandThisLayer(int genomeIndex, int layerIndex, int slotIndex) {

	int relativeExperience =
		genomes[genomeIndex]->objectsToEvaluate -
		ObjectsForLayer(layerIndex);

	while ( relativeExperience>0 ) {
		ExpandLayer(layerIndex,slotIndex);
		RecalculateFrontStrict(layerIndex);
		relativeExperience =
        		genomes[genomeIndex]->objectsToEvaluate -
        		ObjectsForLayer(layerIndex);
	}
}

void   GA::FillSlot(int currentSlot) {

	int genomeToSend = simParams->RandInt(0,simParams->popSize-1);

	FillSlot(genomeToSend,currentSlot);
}

void   GA::FillSlot(int genomeToSend, int currentSlot) {

	int whichObject = int(genomes[genomeToSend]->objectOrdering->Get(0,genomes[genomeToSend]->objectsEvaluated));

	FillSlot(genomeToSend,currentSlot,whichObject);
}

void   GA::FillSlot(int genomeToSend, int currentSlot, int whichObject) {

	GENOME *currentGenome = genomes[genomeToSend];

	SendEarlyStoppingData(currentGenome,currentSlot);
	currentGenome->SendEarlyStoppingData(currentSlot);
	currentGenome->SendAsBrain(-1,currentSlot);
	SendBody(-1,	whichObject,
			currentSlot,currentGenome->objectsToEvaluate);

	genomeIndices[currentSlot] 	= genomeToSend;
	slotStates[currentSlot]    	= RUNNING;

	currentGenome = NULL;
}

int  GA::FindAnyone(int layerIndex) {

	int firstIndexInLayer = simParams->numIndsPerLayer * (layerIndex);
        int lastIndexInLayer  = simParams->numIndsPerLayer * (layerIndex+1) - 1;

        int potentialVictim = simParams->RandInt(firstIndexInLayer,lastIndexInLayer);

	return( potentialVictim );
}

int  GA::FindBest(void) {

	int maxObjects = MostObjectsBeingTackled();

	double maxFitness = -1000.0;
	int bestIndex = -1;

	for (int i=0;i<simParams->popSize;i++)

		if (	genomes[i] &&
			(ObjectsForLayer(AgeLayer(i))==maxObjects) )

			if ( genomes[i]->fitness > maxFitness ) {
				maxFitness = genomes[i]->fitness;
				bestIndex = i;
			}

	return( bestIndex );
}

int  GA::FindBestOnLayer(int layerIndex) {

	double maxFitness = -1000.0;
	int bestIndex = -1;

	for (int i=0;i<simParams->popSize;i++)

		if (	genomes[i] &&
			(AgeLayer(i)==layerIndex) )

			if ( genomes[i]->fitness > maxFitness ) {
				maxFitness = genomes[i]->fitness;
				bestIndex = i;
			}

	return( bestIndex );
}

int  GA::FindDominated(int layerIndex) {

        int firstIndexInLayerAbove = simParams->numIndsPerLayer * (layerIndex);
        int lastIndexInLayerAbove  = simParams->numIndsPerLayer * (layerIndex+1) - 1;
        int foundVictim = false;

        int potentialVictim = simParams->RandInt(firstIndexInLayerAbove,lastIndexInLayerAbove);
        int victimsAssessed = 0;

        while ( (victimsAssessed<simParams->numIndsPerLayer) && (!foundVictim) ) {

                if ( (!Empty(potentialVictim)) && genomes[potentialVictim]->dominated )
                        foundVictim = true;
                else {
                        potentialVictim++;
                        if ( potentialVictim > lastIndexInLayerAbove )
                                potentialVictim = firstIndexInLayerAbove;
                        victimsAssessed++;
                }
        }

        if ( foundVictim )
                return( potentialVictim );
        else
                return( -1 );
}

int GA::FindGenomeOnActivatedLayer(void) {

	int successfulGenome    = FindSuccessfulGenome();
	if ( successfulGenome > -1 )
		return( successfulGenome );

	int validLayer = false;
	int i;

	while ( !validLayer ) {

		    i 			= simParams->RandInt(0,simParams->popSize-1);
		int layerActivated 	= LayerActivated(AgeLayer(i));

		if ( 	layerActivated )
			validLayer = true;
	}

	return( FindAnyone(AgeLayer(i)) );
}

int  GA::FindHole(int layerIndex) {

        int firstIndexInLayerAbove = simParams->numIndsPerLayer * (layerIndex);
        int lastIndexInLayerAbove  = simParams->numIndsPerLayer * (layerIndex+1) - 1;
        int foundVictim = false;

        int potentialVictim = simParams->RandInt(firstIndexInLayerAbove,lastIndexInLayerAbove);
        int victimsAssessed = 0;

        while ( (victimsAssessed<simParams->numIndsPerLayer) && (!foundVictim) ) {

                if ( Empty(potentialVictim) )
                        foundVictim = true;
                else {
                        potentialVictim++;
                        if ( potentialVictim > lastIndexInLayerAbove )
                                potentialVictim = firstIndexInLayerAbove;
                        victimsAssessed++;
                }
        }

        if ( foundVictim )
                return( potentialVictim );
        else
                return( -1 );
}

int GA::FindInferior(int aggressor, int layerIndex) {

	int found = false;
	int victimIndex=0;

	while( (!found) && (victimIndex<simParams->popSize) ) {

		if ( 	(genomes[victimIndex]) &&
			(AgeLayer(victimIndex)==layerIndex) &&
			(genomes[aggressor]->SuperiorTo(genomes[victimIndex]) ) ) {

			found = true;
		}
		else {
			victimIndex++;
		}
	}

	if ( found )
		return( victimIndex );
	else
		return( -1 );
}

int  GA::FindLayerForNewChild(void) {

/*
	int successfulGenome = FindSuccessfulGenome();
	if ( successfulGenome>-1 )
		return( AgeLayer(successfulGenome) );
*/
	return( simParams->RandInt(0,simParams->numLayers-1) );
}

int    GA::FindNonDominatedAndNotTooOld(int layerIndex) {

        int firstIndexInLayerAbove = simParams->numIndsPerLayer * (layerIndex);
        int lastIndexInLayerAbove  = simParams->numIndsPerLayer * (layerIndex+1) - 1;
        int foundVictim = false;

        int potentialVictim = simParams->RandInt(firstIndexInLayerAbove,lastIndexInLayerAbove);
        int victimsAssessed = 0;

        while ( (victimsAssessed<simParams->numIndsPerLayer) && (!foundVictim) ) {

                if ( 	(!Empty(potentialVictim)) && 
			(!genomes[potentialVictim]->dominated) &&
			(!TooOld(potentialVictim)) )
                        foundVictim = true;
                else {
                        potentialVictim++;
                        if ( potentialVictim > lastIndexInLayerAbove )
                                potentialVictim = firstIndexInLayerAbove;
                        victimsAssessed++;
                }
        }

        if ( foundVictim )
                return( potentialVictim );
        else
                return( -1 );
}

int    GA::FindNonDominatedAndTooOld(int layerIndex) {

        int firstIndexInLayerAbove = simParams->numIndsPerLayer * (layerIndex);
        int lastIndexInLayerAbove  = simParams->numIndsPerLayer * (layerIndex+1) - 1;
        int foundVictim = false;

        int potentialVictim = simParams->RandInt(firstIndexInLayerAbove,lastIndexInLayerAbove);
        int victimsAssessed = 0;

        while ( (victimsAssessed<simParams->numIndsPerLayer) && (!foundVictim) ) {

                if (    (!Empty(potentialVictim)) &&
                        (!genomes[potentialVictim]->dominated) &&
                        (TooOld(potentialVictim)) )
                        foundVictim = true;
                else {
                        potentialVictim++;
                        if ( potentialVictim > lastIndexInLayerAbove )
                                potentialVictim = firstIndexInLayerAbove;
                        victimsAssessed++;
                }
        }

        if ( foundVictim )
                return( potentialVictim );
        else
                return( -1 );
}

int    GA::FindParent(int layerIndex) {

	return( FindNonDominatedAndNotTooOld(layerIndex) );
/*
        int potentialParent = FindNonDominatedAndNotTooOld(layerIndex);
        if ( potentialParent != -1 )
                return( potentialParent );

        potentialParent = FindNonDominatedAndTooOld(layerIndex);

        return( potentialParent );
*/

}

int    GA::FindSlot(void) {

	return( simParams->RandInt(0,simParams->numSlots-1) );
}

int GA::FindSuccessful(int layerIndex) {

        int firstIndexInLayerAbove = simParams->numIndsPerLayer * (layerIndex);
        int lastIndexInLayerAbove  = simParams->numIndsPerLayer * (layerIndex+1) - 1;
        int foundVictim = false;

        int potentialVictim = simParams->RandInt(firstIndexInLayerAbove,lastIndexInLayerAbove);
        int victimsAssessed = 0;

        while ( (victimsAssessed<simParams->numIndsPerLayer) && (!foundVictim) ) {

                if ( (!Empty(potentialVictim)) && genomes[potentialVictim]->Successful() )
                        foundVictim = true;
                else {
                        potentialVictim++;
                        if ( potentialVictim > lastIndexInLayerAbove )
                                potentialVictim = firstIndexInLayerAbove;
                        victimsAssessed++;
                }
        }

        if ( foundVictim )
                return( potentialVictim );
        else
                return( -1 );
}

int GA::FindSuccessfulGenome(void) {

	int successfulGenome = -1;

	int i=0;

	while ( (i<simParams->popSize) && (successfulGenome==-1) )
		if ( genomes[i] && genomes[i]->Successful() )
			successfulGenome = i;
		else
			i++;

	return( successfulGenome );
}

int GA::FindTarget(int layerIndex) {

	int potentialTarget = FindTooOld(layerIndex);
	if ( potentialTarget != -1 )
		return( potentialTarget );

	potentialTarget = FindHole(layerIndex);
	if ( potentialTarget != -1 )
		return potentialTarget;

	potentialTarget = FindDominated(layerIndex);
	return( potentialTarget ); 
}

int GA::FindTooOld(int layerIndex) {

        int firstIndexInLayerAbove = simParams->numIndsPerLayer * (layerIndex);
        int lastIndexInLayerAbove  = simParams->numIndsPerLayer * (layerIndex+1) - 1;
        int foundVictim = false;

        int potentialVictim = simParams->RandInt(firstIndexInLayerAbove,lastIndexInLayerAbove);
        int victimsAssessed = 0;

        while ( (victimsAssessed<simParams->numIndsPerLayer) && (!foundVictim) ) {

                if ( (!Empty(potentialVictim)) && TooOld(potentialVictim) )
                        foundVictim = true;
                else {
                        potentialVictim++;
                        if ( potentialVictim > lastIndexInLayerAbove )
                                potentialVictim = firstIndexInLayerAbove;
                        victimsAssessed++;
                }
        }

        if ( foundVictim )
                return( potentialVictim );
        else
                return( -1 );
}

int GA::FindTooOldGenome(void) {

	int found = false;
	int potentialLayer = simParams->RandInt(0,simParams->numLayers-1);
	int layersAssessed = 0;

	while (	(!found) &&
		(layersAssessed<simParams->numLayers) ) {

		int tooOld = FindTooOld(potentialLayer);
		if ( tooOld > -1 )
			found = true;
		else {
			potentialLayer++;
			if ( potentialLayer == simParams->numLayers )
				potentialLayer = 0;
			layersAssessed++;
		}
	}

	if ( found )
		return( potentialLayer );
	else
		return( -1 );
}

int GA::FindVictim(int i, int layerIndex) {

        int potentialVictim = FindSuccessful(layerIndex);
        if ( potentialVictim != -1 )
                return( potentialVictim );

        potentialVictim = FindTooOld(layerIndex);
        if ( potentialVictim != -1 )
                return( potentialVictim );

        potentialVictim = FindHole(layerIndex);
        if ( potentialVictim != -1 )
                return( potentialVictim );

	return( FindInferior(i,layerIndex) );
	
/*
        potentialVictim = FindInferior(i,layerIndex);
        if ( potentialVictim != -1 )
                return( potentialVictim );

        int tooExperienced = TooMuchExperienceForLayer(i,layerIndex);
        int nonDominated   = (!Dominated(genomes[i],layerIndex));
        int equalExp       = genomes[i]->objectsEvaluated == ObjectsForLayer(layerIndex);

        int strongerThanDominated =     tooExperienced ||
                                        (nonDominated && equalExp);

        if ( strongerThanDominated )
                return( FindDominated(layerIndex) );
        else
                return( -1 );
*/

}

void GA::FixLayer(int layerIndex, int slotIndex) {

	RecalculateFrontStrict(layerIndex);

	while ( NumNonDominated(layerIndex)==0 ) {

		int mostPromisingGenome = BestGenomeOnLayer(layerIndex);

		EvaluateGenomeOnce(mostPromisingGenome,slotIndex);
		RecalculateFrontStrict(layerIndex);
	}
}

void GA::FixLayers(int slotIndex) {

	for (int i=0;i<simParams->numLayers;i++)
		FixLayer(i,slotIndex);
}

int  GA::GenomesToEvaluateOnThisLayer(int layerIndex) {
	
	if ( NumNonDominated(layerIndex)==0 )
		return( true );

        int found = false;
        int i=0;

        while ( (!found) &&
                (i<simParams->popSize) ) {

                if (    (genomes[i]) &&
                        (AgeLayer(i)==layerIndex) &&
                        (genomes[i]->objectsEvaluated < genomes[i]->objectsToEvaluate) &&
                        (!Dominated(genomes[i],layerIndex)) ) {

                        found = true;
                }
                else {
                        i++;
                }
        }

        return( found );
}

double GA::GetFingerAngle(int fingerIndex, MATRIX *fingerSpacings) {

	return( fingerSpacings->Get(0,fingerIndex) );
}

int GA::GetSlotState(int currentSlot) {

	if ( slotStates[currentSlot]==EMPTY )
		return( EMPTY );

	else if ( slotStates[currentSlot]==FINISHED )
		return( FINISHED );

	else if ( slotStates[currentSlot]==RUNNING ) {

		char fileName[200];

		sprintf(fileName,"/tmp/Files%d_%d/myCoreDists.dat",simParams->randSeed,currentSlot);
		if ( !simParams->FileExists(fileName) )
			return( RUNNING );

		sprintf(fileName,"/tmp/Files%d_%d/myLeftShoulderDists.dat",simParams->randSeed,currentSlot);
		if ( !simParams->FileExists(fileName) )
			return( RUNNING );

		sprintf(fileName,"/tmp/Files%d_%d/myRightShoulderDists.dat",simParams->randSeed,currentSlot);
		if ( !simParams->FileExists(fileName) )
			return( RUNNING );
		
		sprintf(fileName,"/tmp/Files%d_%d/Sensors.dat",simParams->randSeed,currentSlot);
		if ( !simParams->FileExists(fileName) )
			return( RUNNING );

		slotStates[currentSlot] = FINISHED;
		return( FINISHED );
	}
}

int  GA::HandleAggressorDominatingVictim(int aggressor, int victimIndex, int slotIndex) {

	if ( 	genomes[aggressor]->objectsEvaluated >= 
		genomes[victimIndex]->objectsEvaluated )

		return( true );	// Victim failed

	while ( genomes[aggressor]->objectsToEvaluate <
		genomes[victimIndex]->objectsToEvaluate )
		genomes[aggressor]->Expand();

	while ( genomes[aggressor]->Dominates(genomes[victimIndex]) &&
		(genomes[aggressor]->objectsEvaluated < 
		 genomes[victimIndex]->objectsEvaluated) ) {

		EvaluateGenomeOnce(aggressor,slotIndex);
		char note[100];
		sprintf(note,"Expanded genome.");
		PrintGenome(aggressor,note);
	}

	if ( genomes[aggressor]->Dominates(genomes[victimIndex]) )
		return( true );	// Victim failed

	return( false );	// Aggressor failed
}

int  GA::HandleNeitherDominatingOneAnother(int aggressor, int victimIndex) {

	int aggressorHasMoreExperience = 	genomes[aggressor]->objectsEvaluated > 
						genomes[victimIndex]->objectsEvaluated;

	int victimFailed = aggressorHasMoreExperience;

	return( victimFailed );
}

void GA::HandleSuccessfulGenome(int successfulIndex, int slotIndex) {

	if ( (successfulIndex>-1) && genomes[successfulIndex]->Successful() ) {

//		if ( genomes[successfulIndex]->objectsEvaluated==simParams->totalObjects ) {
		if ( 	(simParams->noScaffolding) ||
			(simParams->stanceIndex == (simParams->totalObjects-1)) ) {
			WriteAll();
			exit(0);
		}

		SaveBest();

		simParams->successOccurred = false;
		simParams->successfulGenomeIndex = -1;

		EmptyAllSlots();

		simParams->stanceIndex++;

		SortGenomesOnLayers();

		ResetAllGenomes();

		EvaluateAllLayers();
	}
}

int  GA::InferiorOnLayer(int genomeIndex, int layerIndex) {

	int inferior = false;

	int firstIndexInLayer = simParams->numIndsPerLayer * (layerIndex);
        int lastIndexInLayer  = simParams->numIndsPerLayer * (layerIndex+1) - 1;

	int otherIndex = firstIndexInLayer;

	while ( (otherIndex<=lastIndexInLayer) && (!inferior) ) {

		if ( 	genomes[otherIndex] &&
			genomes[otherIndex]->SuperiorTo(genomes[genomeIndex]) )
			inferior = true;
		else
			otherIndex++;
	}

	return( inferior );
}

int GA::InTopLayer(int genomeIndex) {

	return( AgeLayer(genomeIndex)==(simParams->numLayers-1) );
}

int  GA::InvalidGenomeOnLayer(int layerIndex) {

	int invalidGenomeIndex = -1;

	for (int i=0;i<simParams->popSize;i++)

		if ( 	(AgeLayer(i)==layerIndex) &&
			(genomes[i]) &&
			(!genomes[i]->dominated) &&
			(genomes[i]->objectsEvaluated<(layerIndex+1)) )

			invalidGenomeIndex = i;

	return( invalidGenomeIndex );
}

int GA::IsAVictim(int potentialVictim, int aggressor) {

	if (	Empty(potentialVictim) )
		return( true );

	int	tooOld	  	= TooOld(potentialVictim);

	int	inferior 	= genomes[aggressor]->SuperiorTo(genomes[potentialVictim]);

//	if ( 	InTopLayer(potentialVictim) || 
//		(!LayerActivated(AgeLayer(potentialVictim)+1)) )
//		return( inferior );
//	else
		return( tooOld || inferior );
}

int  GA::IsNonDominatedInLayer(int genomeIndex, int ageLayer) {

	int isNonDominated = true;

	int otherGenomeIndex=0;

	while ( 	(otherGenomeIndex<simParams->popSize) && 
			(isNonDominated) ) {

		if ( 	(genomeIndex!=otherGenomeIndex) &&
			(AgeLayer(otherGenomeIndex)==ageLayer) )

			isNonDominated = genomes[otherGenomeIndex]->Dominates(genomes[genomeIndex]);

		otherGenomeIndex++;
	}

	return( isNonDominated );
}

void GA::KillDominated(void) {

	for (int i=0;i<simParams->popSize;i++)

		if ( genomes[i]->dominated ) {

			SeverChildConnections(genomes[i]);
			delete genomes[i];
			genomes[i] = NULL;
		}
}

int  GA::LayerActivated(int layerIndex) {

	int activated = false;

	for (int i=0;i<simParams->popSize;i++)
		if ( AgeLayer(i)==layerIndex )
			if ( genomes[i] )
				activated = true;

	return( activated );
}

int  GA::LayerEvaluated(int layerIndex) {

	int evaluated = true;

	int i=0;

	while ( (evaluated) && (i<simParams->popSize) ) {

		if ( AgeLayer(i)==layerIndex ) {
			if ( genomes[i]->objectsEvaluated==0 )
				evaluated = false;
			else
				i++;
		}
		else
			i++;
	}

	return( evaluated );
}

int    GA::MaxAgeForLayer(int layerIndex) {

	return( AGE_CEILINGS[layerIndex] );
}

int GA::MaxEvalsOnLayer(int layerIndex) {

        int maxEvals = -1;

        for (int i=0;i<simParams->popSize;i++)

                if (    (genomes[i]) &&
                        (AgeLayer(i)==layerIndex) &&
                        (genomes[i]->objectsEvaluated > maxEvals) )

                        maxEvals = genomes[i]->objectsEvaluated;

        return( maxEvals );
}

double GA::MaxMult(void) {

	double maxMult = -1;
	int    maxMultIndex = 0;

	for (int i=0;i<simParams->popSize;i++) {
		if ( genomes[i] )
			if ( (genomes[i]->coreDist * genomes[i]->leftShoulderDist * genomes[i]->rightShoulderDist) > maxMult ) {
				maxMult = genomes[i]->coreDist * genomes[i]->leftShoulderDist * genomes[i]->rightShoulderDist;
				maxMultIndex = i;
			}
	}

	return( maxMult );
}

int  GA::MaxMultIndex(void) {

	double maxMult = -1;
	int    maxMultIndex = 0;

	for (int i=0;i<simParams->popSize;i++) {
		if ( genomes[i] )
			if ( (genomes[i]->coreDist * genomes[i]->leftShoulderDist * genomes[i]->rightShoulderDist) > maxMult ) {
				maxMult = genomes[i]->coreDist * genomes[i]->leftShoulderDist * genomes[i]->rightShoulderDist;
				maxMultIndex = i;
			}
	}

	return( maxMultIndex );
}

int  GA::MaxSuccessfulObjectsOnLayer(int layerIndex) {

	int maxSuccessfulObjects = 0;

	for (int i=0;i<simParams->popSize;i++) {

		if (	(genomes[i])	&&
			(AgeLayer(i)==layerIndex) ) {

			int mo = 0;

			for (int j=0;j<genomes[i]->objectsEvaluated;j++) {

				if ( genomes[i]->SuccessfulOnObject(j) )
					mo++;
			}
		
/*	
			while (	(genomes[i]->SuccessfulOnObject(mo))	&&
				(genomes[i]->objectsEvaluated>mo) )
				mo++;
*/
			if ( mo > maxSuccessfulObjects )

				maxSuccessfulObjects = mo;
		}
	}

	return( maxSuccessfulObjects );
}

int GA::MostObjectsBeingTackled(void) {

	int mostObjects = -1;

	for (int i=0;i<simParams->numLayers;i++)

		if ( ObjectsForLayer(i) > mostObjects )
			mostObjects = ObjectsForLayer(i);

	return( mostObjects );
}

void GA::Move(int from, int to) {

	if ( genomes[from]->Successful() ) {
		simParams->successfulGenomeIndex = to;
	}

	genomes[to] = genomes[from];
	genomes[from] = NULL;

	char note[100];
	sprintf(note,"Moved genome from %d to %d.",from,to);
	PrintGenome(to,note);	

//	if ( to != simParams->popSize )
//		RecalculateFrontStrict(AgeLayer(to));

	
}

int GA::MovedLongEnoughAgo(int genomeIndex) {

	if ( !genomes[genomeIndex] )
		return( true );

	int timeSinceLastMove = int(numGenomesEvaluated) - int(genomes[genomeIndex]->evalsAtMoveUp); 

	return( timeSinceLastMove > simParams->popSize );
}

void GA::MoveReinitializationPointer(void) {

	reinitializationIndividual++;
}

int  GA::NoSuccessesForObject(int objectIndex, int layerIndex) {

        int found = false;
        int i = 0;

        while ( (!found) && (i<simParams->popSize) ) {

                if (    (genomes[i]) &&
                        (AgeLayer(i)==layerIndex) &&
                        (genomes[i]->objectsEvaluated>=(objectIndex+1)) &&
                        (genomes[i]->SuccessfulOnObject(objectIndex)) )

                        found = true;
                else
                        i++;
        }

        int noSuccessesFound = (found==false);

        return( noSuccessesFound );
}

int  GA::NotEnoughExperienceForLayer(int genomeIndex, int layerIndex) {

	return( genomes[genomeIndex]->objectsToEvaluate < ObjectsForLayer(layerIndex) );
}

int  GA::NumFilled(int ageLayer) {

	int num = 0;

	for (int i=0;i<simParams->popSize;i++) {
		if ( genomes[i] && (AgeLayer(i)==ageLayer) )
			num++;
	}

	return( num );
}

int  GA::NumNonDominated(int ageLayer) {

	int num = 0;

	for (int i=0;i<simParams->popSize;i++)
		if ( 	(genomes[i]) && 
			(AgeLayer(i)==ageLayer) && 
			(!genomes[i]->dominated) )
			num++;

	return( num );
}

int  GA::NumNonDominated(GENOME *currGenome, int ageLayer) {

	int num = 0;

	for (int i=0;i<simParams->popSize;i++)
		if ( 	(genomes[i]) && 
			(AgeLayer(i)==ageLayer) && 
			(!genomes[i]->dominated) && 
			(genomes[i]!=currGenome) )
			num++;

	return( num );
}

int  GA::ObjectsForLayer(int layerIndex) {

	return( objectsForLayer[layerIndex] );
}

void GA::PrintGenome(int genomeIndex, char *note) {

	printf("[nge: %0.0f] ",numGenomesEvaluated);

	printf("[gi: %d] ",genomeIndex);
	if ( genomeIndex<100 )
		printf(" ");
	if ( genomeIndex<10 )
		printf(" ");

	printf("[age: %d %d %d] ",
		AgeLayer(genomeIndex),
		genomes[genomeIndex]->Age(int(numGenomesEvaluated)),
		MaxAgeForLayer(AgeLayer(genomeIndex)));

	if ( AgeLayer(genomeIndex)<10 )
		printf(" ");

	if ( genomes[genomeIndex]->Age(int(numGenomesEvaluated))<10 )
		printf(" ");

        if ( MaxAgeForLayer(AgeLayer(genomeIndex))<10000 )
                printf(" ");
        if ( MaxAgeForLayer(AgeLayer(genomeIndex))<1000 )
                printf(" ");
        if ( MaxAgeForLayer(AgeLayer(genomeIndex))<100 )
                printf(" ");
	if ( MaxAgeForLayer(AgeLayer(genomeIndex))<10 )
		printf(" ");

	printf("[l:");
	for (int i=0;i<simParams->numLayers;i++) {
//		int objectsForLayer = ObjectsForLayer(i);
//		printf("%d",objectsForLayer);
		printf("%d",simParams->stanceIndex);
	}
	printf("] ");

	printf("[dom: %d %d] ",
		genomes[genomeIndex]->dominated,
		NumNonDominated(AgeLayer(genomeIndex)));

	printf("[obj: %d %d %d] ",
		ObjectsForLayer(AgeLayer(genomeIndex)),
		genomes[genomeIndex]->objectsEvaluated,
		genomes[genomeIndex]->objectsToEvaluate);

	printf("[fit: %5.5f] ",
		genomes[genomeIndex]->fitness);
//        printf("[fit: %3.3f %3.3f %3.3f] ",genomes[genomeIndex]->coreDist,
  //                                      genomes[genomeIndex]->leftShoulderDist,
    //                                    genomes[genomeIndex]->rightShoulderDist);

//	int bestIndex = FindBestOnLayer(AgeLayer(genomeIndex));
	int bestIndex = FindBest();

	printf("[best: %d %d %d %d ",
		bestIndex,
		AgeLayer(bestIndex),
		genomes[bestIndex]->Age(int(numGenomesEvaluated)),
		ObjectsForLayer(AgeLayer(bestIndex)));

	printf("%5.5f] ",
		genomes[bestIndex]->fitness);
//        printf("%3.3f %3.3f %3.3f] ",
//                genomes[bestIndex]->coreDist,
//                genomes[bestIndex]->leftShoulderDist,
//                genomes[bestIndex]->rightShoulderDist);

	printf("%s \n",note);
}

void GA::PrintLayer(int layerIndex) {

	char note[100];
	sprintf(note,"Checking genome.");

	for (int i=0;i<simParams->popSize;i++)

		if ( 	genomes[i] &&
			(AgeLayer(i)==layerIndex) )

			PrintGenome(i,note);
}

void GA::PrintSlots(void) {

	for (int i=0;i<simParams->numSlots;i++)
		if ( slotStates[i]==EMPTY )
			printf("E \t");
		else if ( slotStates[i]==RUNNING )
			printf("R \t");
		else if ( slotStates[i]==FINISHED )
			printf("F \t");
	printf("\n");

	for (int i=0;i<simParams->numSlots;i++)
		printf("%d \t",genomeIndices[i]);
	printf("\n");

}

int  GA::RandomGenome(int ageLayer) {

	int firstIndexOfLayer = simParams->numIndsPerLayer * ageLayer;
	int lastIndexOfLayer  = simParams->numIndsPerLayer * (ageLayer+1) - 1;

	return( simParams->RandInt(firstIndexOfLayer,lastIndexOfLayer) );
}

void GA::RandomizeGenome(int genomeIndex) {

	if ( genomes[genomeIndex] ) {
		delete genomes[genomeIndex];
		genomes[genomeIndex] = new GENOME(int(numGenomesEvaluated));
	}
}

void   GA::RecalculateFront(int ageLayer) {

	for (int i=0;i<simParams->popSize;i++) {

		if ( 	genomes[i] && 
			(AgeLayer(i)==ageLayer) ) 

			genomes[i]->dominated = false;
	}

	for (int i=0;i<simParams->popSize;i++) {

		int j=0;
		while ( 	(j<simParams->popSize) && 
				(genomes[i]) &&
				(AgeLayer(i)==ageLayer) &&
				(!genomes[i]->dominated) ) {

			if ( 	(i!=j) && 
				(genomes[j]) &&
				(AgeLayer(j)==ageLayer) )

				genomes[i]->dominated = genomes[j]->SuperiorTo(genomes[i]);

			j++;
		}
	}
}

void   GA::RecalculateFrontStrict(int ageLayer) {

	for (int i=0;i<simParams->popSize;i++) {

		if ( 	genomes[i] && 
			(AgeLayer(i)==ageLayer) ) 

			if ( genomes[i]->objectsEvaluated<ObjectsForLayer(ageLayer) )
				genomes[i]->dominated = true;
			else
				genomes[i]->dominated = false;
	}

	for (int i=0;i<simParams->popSize;i++) {

		int j=0;
		while ( 	(j<simParams->popSize) && 
				(genomes[i]) &&
				(AgeLayer(i)==ageLayer) &&
				(!genomes[i]->dominated) ) {

			if ( 	(i!=j) && 
				(genomes[j]) &&
				(AgeLayer(j)==ageLayer) )

				genomes[i]->dominated = genomes[j]->Superior(genomes[i]);

			j++;
		}
	}
}

void   GA::RecalculateFrontsStrict(int victimLayer) {

	for (int i=victimLayer;i<simParams->numLayers;i++)

		RecalculateFrontStrict(i);
}

void   GA::RecalculateFrontUsingDomination(int ageLayer) {

	for (int i=0;i<simParams->popSize;i++) {

		if ( 	genomes[i] && 
			(AgeLayer(i)==ageLayer) ) 

			genomes[i]->dominated = false;
	}

	for (int i=0;i<simParams->popSize;i++) {

		int j=0;
		while ( 	(j<simParams->popSize) && 
				(genomes[i]) &&
				(AgeLayer(i)==ageLayer) &&
				(!genomes[i]->dominated) ) {

			if ( 	(i!=j) && 
				(genomes[j]) &&
				(AgeLayer(j)==ageLayer) )

				genomes[i]->dominated = genomes[j]->Dominates(genomes[i]);

			j++;
		}
	}
}

int GA::ReinitializationFinished(void) {

	return( reinitializationIndividual==simParams->numIndsPerLayer );
}

void GA::ReplaceVictimAbove(int i, int slotIndex) {

        int victimIndex = FindVictim(i,AgeLayer(i)+1);
	if ( victimIndex==-1 ) {
		Discard(i);
		return;
	}

        TryMoveUp(victimIndex,slotIndex);

	if ( !Empty(victimIndex) )
		
		Discard(victimIndex);

	Move(i,victimIndex);
/*
	int victimIndex = FindVictim(i,AgeLayer(i)+1);

	if ( victimIndex==-1 ) {
		if (    TooMuchExperienceForLayer(i,AgeLayer(i)) )
			Discard(i);
		return;
	}

	TryMoveUp(victimIndex,slotIndex);

	printf("Victim failed?\n");
	int victimFailed 	= VictimFailed(i, victimIndex,slotIndex);
	printf("Victim failure: %d\n",victimFailed);

	if ( victimFailed ) {

		if ( !Empty(victimIndex) )
			Discard(victimIndex);

		printf("About to move %d.\n",i);
        	Move(i,victimIndex);
		printf("Moved %d.\n",i);

		printf("About to balance layer.\n");
		BalanceLayer(victimIndex,AgeLayer(victimIndex),slotIndex);
		printf("Balanced layer.\n");

        	char note[100];
        	sprintf(note,"genome moved from %d to %d.",i,victimIndex);
        	PrintGenome(victimIndex,note);
	}
	else {
		if (	TooMuchExperienceForLayer(i,AgeLayer(i)) )
			Discard(i);
	}
*/
}

void GA::ResetAllGenomes(void) {

	for (int i=0;i<simParams->popSize;i++)

		genomes[i]->Reset();
}

void GA::SaveBest(void) {

	int bestIndex = FindBest();

	int whichObject = int(genomes[bestIndex]->objectOrdering->Get(0,0));

	genomes[bestIndex]->SendAsBrain(simParams->stanceIndex,-1);
	SendBody(simParams->stanceIndex,	whichObject,-1,genomes[bestIndex]->objectsToEvaluate);

	previousBestFit = currBestFit;

	printf("%3.3f %3.3f %3.3f\n",	genomes[bestIndex]->coreDist,
					genomes[bestIndex]->leftShoulderDist,
					genomes[bestIndex]->rightShoulderDist);

	currBestFit = 			genomes[bestIndex]->coreDist * 
					genomes[bestIndex]->leftShoulderDist * 
					genomes[bestIndex]->rightShoulderDist;
}

void GA::SaveObjsConquered(void) {

	MATRIX *objsConq = new MATRIX(1,simParams->popSize,0.0);

	for (int j=0;j<simParams->popSize;j++)
		objsConq->Set(0,j,genomes[j]->objectsEvaluated);

	objsConq->AppendToFile(simParams->objsConqFileName);

	delete objsConq;
	objsConq = NULL;
}

void GA::SaveOutData(void) {

	if ( rowIndex == (DATA_FILE_BUFFER-1) ) {
		dataForOutputFile->AppendToFile(dataFileName);
		dataForOutputFile->ReZero();
		rowIndex=0;
	}
	else
		rowIndex++;
}

void GA::SendArm_Joints(ofstream *outFile) {

	(*outFile) << "\n( FrontSegment_LeftShoulderGirdle\n";
	(*outFile) << "-connect FrontSegment LeftShoulderGirdle\n";
	(*outFile) << "-jointPosition "	<< -SEGMENT_WIDTH/2.0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< +SEGMENT_LENGTH - 2.0*CORE_RADIUS<<" \n";
	(*outFile) << "-jointType Prismatic\n";
	(*outFile) << "-jointNormal -1 0 0\n";
	(*outFile) << "-jointLimits "	<< startRetraction << " "
					<< endExtension << " "
					<< endStanceTransitionWhen << "\n";
	(*outFile) << ")\n";

	(*outFile) << "\n( FrontSegment_RightShoulderGirdle\n";
	(*outFile) << "-connect FrontSegment RightShoulderGirdle\n";
	(*outFile) << "-jointPosition "	<< SEGMENT_WIDTH/2.0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< +SEGMENT_LENGTH - 2.0*CORE_RADIUS<<" \n";
	(*outFile) << "-jointType Prismatic\n";
	(*outFile) << "-jointNormal 1 0 0\n";
	(*outFile) << "-jointLimits "	<< startRetraction << " "
					<< endExtension << " "
					<< endStanceTransitionWhen << "\n";
	(*outFile) << ")\n";



	(*outFile) << "\n( LeftShoulderGirdle_LeftShoulder\n";
	(*outFile) << "-connect LeftShoulderGirdle LeftShoulder\n";
	(*outFile) << "-jointPosition "	<< -SEGMENT_WIDTH/2.0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< SEGMENT_LENGTH - 2.0*CORE_RADIUS<<" \n";
	(*outFile) << "-jointType Hinge\n";
	(*outFile) << "-jointNormal 0 0 1\n";
	(*outFile) << "-jointLimits "	<<	stance->Get(0,0)*ratio
					<<" "<<	stance->Get(1,0)*ratio
					<<" "<< endStance->Get(0,0)*endRatio
					<<" "<<	endStance->Get(1,0)*endRatio
					<<" "<< endStanceTransitionWhen<<" \n";
	(*outFile) << "-addMotor\n";
	(*outFile) << "-motorForce "<<MOTOR_FORCE<<"\n";
	(*outFile) << "-motorSpeed "<<MOTOR_SPEED<<"\n";
        (*outFile) << "-addSensor \n";
	(*outFile) << ")\n";

	(*outFile) << "\n( RightShoulderGirdle_RightShoulder\n";
	(*outFile) << "-connect RightShoulderGirdle RightShoulder\n";
	(*outFile) << "-jointPosition "	<< SEGMENT_WIDTH/2.0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< SEGMENT_LENGTH - 2.0*CORE_RADIUS<<" \n";
	(*outFile) << "-jointType Hinge\n";
	(*outFile) << "-jointNormal 0 0 -1\n";
	(*outFile) << "-jointLimits "	<<	stance->Get(0,0)*ratio
					<<" "<<	stance->Get(1,0)*ratio
					<<" "<< endStance->Get(0,0)*endRatio
					<<" "<<	endStance->Get(1,0)*endRatio
					<<" "<< endStanceTransitionWhen<<" \n";
	(*outFile) << "-addMotor\n";
	(*outFile) << "-motorForce "<<MOTOR_FORCE<<"\n";
	(*outFile) << "-motorSpeed "<<MOTOR_SPEED<<"\n";
        (*outFile) << "-addSensor \n";
	(*outFile) << ")\n";


	(*outFile) << "\n( LeftShoulder_LeftArm\n";
	(*outFile) << "-connect LeftShoulder LeftArm\n";
	(*outFile) << "-jointPosition "	<< -SEGMENT_WIDTH/2.0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< SEGMENT_LENGTH - 2.0*CORE_RADIUS<<" \n";
	(*outFile) << "-jointType Hinge\n";
	(*outFile) << "-jointNormal 0 1 0\n";
	(*outFile) << "-jointLimits "	<<	stance->Get(0,1)*ratio
					<<" "<<	stance->Get(1,1)*ratio
					<<" "<< endStance->Get(0,1)*endRatio
					<<" "<<	endStance->Get(1,1)*endRatio
					<<" "<< endStanceTransitionWhen<<" \n";
	(*outFile) << "-addMotor\n";
	(*outFile) << "-motorForce "<<MOTOR_FORCE<<"\n";
	(*outFile) << "-motorSpeed "<<MOTOR_SPEED<<"\n";
        (*outFile) << "-addSensor \n";
	(*outFile) << ")\n";

	(*outFile) << "\n( RightShoulder_RightArm\n";
	(*outFile) << "-connect RightShoulder RightArm\n";
	(*outFile) << "-jointPosition "	<< SEGMENT_WIDTH/2.0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< SEGMENT_LENGTH - 2.0*CORE_RADIUS<<" \n";
	(*outFile) << "-jointType Hinge\n";
	(*outFile) << "-jointNormal 0 -1 0\n";
	(*outFile) << "-jointLimits "	<<	stance->Get(0,1)*ratio
					<<" "<<	stance->Get(1,1)*ratio
					<<" "<< endStance->Get(0,1)*endRatio
					<<" "<<	endStance->Get(1,1)*endRatio
					<<" "<< endStanceTransitionWhen<<" \n";
	(*outFile) << "-addMotor\n";
	(*outFile) << "-motorForce "<<MOTOR_FORCE<<"\n";
	(*outFile) << "-motorSpeed "<<MOTOR_SPEED<<"\n";
        (*outFile) << "-addSensor \n";
	(*outFile) << ")\n";


	(*outFile) << "\n( LeftArm_LeftLowerArm\n";
	(*outFile) << "-connect LeftArm LeftLowerArm\n";
	(*outFile) << "-jointPosition "	<< -SEGMENT_WIDTH/2.0 -ARM_LENGTH/2.0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< +SEGMENT_LENGTH - 2.0*CORE_RADIUS<<" \n";
	(*outFile) << "-jointType Prismatic\n";
	(*outFile) << "-jointNormal -1 0 0\n";
	(*outFile) << "-jointLimits "	<< startRetraction << " "
					<< endExtension << " "
					<< endStanceTransitionWhen << "\n";
	(*outFile) << ")\n";

	(*outFile) << "\n( RightArm_RightLowerArm\n";
	(*outFile) << "-connect RightArm RightLowerArm\n";
	(*outFile) << "-jointPosition "	<< SEGMENT_WIDTH/2.0 +ARM_LENGTH/2.0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< +SEGMENT_LENGTH - 2.0*CORE_RADIUS<<" \n";
	(*outFile) << "-jointType Prismatic\n";
	(*outFile) << "-jointNormal 1 0 0\n";
	(*outFile) << "-jointLimits "	<< startRetraction << " "
					<< endExtension << " "
					<< endStanceTransitionWhen << "\n";
	(*outFile) << ")\n";
}

void GA::SendArms(ofstream *outFile) {

	(*outFile) << "\n( LeftArm \n";
	(*outFile) << "-position " << -SEGMENT_WIDTH/2.0-ARM_LENGTH/4.0 << " "
				   << SEGMENT_HEIGHT/2.0 << " "
				   << +SEGMENT_LENGTH - 2.0*CORE_RADIUS<< " "
				   << "\n";
	(*outFile) << "-shape -cylinder "
		   << SEGMENT_HEIGHT/2.0 <<" "
		   << ARM_LENGTH/2.0 <<"\n";
	(*outFile) << "-rotation 1 0 0\n";
	(*outFile) << "-mass 0.5 \n";
	(*outFile) << "-colour 0 1 0 \n";
	(*outFile) << ") \n";

	(*outFile) << "\n( RightArm \n";
	(*outFile) << "-position " << +SEGMENT_WIDTH/2.0+ARM_LENGTH/4.0 << " "
				   << SEGMENT_HEIGHT/2.0 << " "
				   << +SEGMENT_LENGTH - 2.0*CORE_RADIUS<< " "
				   << "\n";
	(*outFile) << "-shape -cylinder "
		   << SEGMENT_HEIGHT/2.0 <<" "
		   << ARM_LENGTH/2.0 <<"\n";
	(*outFile) << "-rotation 1 0 0\n";
	(*outFile) << "-mass 0.5 \n";
	(*outFile) << "-colour 0 0.5 0 \n";
	(*outFile) << ") \n";


	(*outFile) << "\n( LeftLowerArm \n";
	(*outFile) << "-position " << -SEGMENT_WIDTH/2.0-3.0*ARM_LENGTH/4.0 << " "
				   << SEGMENT_HEIGHT/2.0 << " "
				   << +SEGMENT_LENGTH - 2.0*CORE_RADIUS<< " "
				   << "\n";
	(*outFile) << "-shape -cylinder "
		   << SEGMENT_HEIGHT/2.0 <<" "
		   << ARM_LENGTH/2.0 <<"\n";
	(*outFile) << "-rotation 1 0 0\n";
	(*outFile) << "-mass 0.5 \n";
	(*outFile) << "-colour 0 1 0 \n";
	(*outFile) << "-addTouchSensor \n";
	(*outFile) << ") \n";

	(*outFile) << "\n( RightLowerArm \n";
	(*outFile) << "-position " << +SEGMENT_WIDTH/2.0+3.0*ARM_LENGTH/4.0 << " "
				   << SEGMENT_HEIGHT/2.0 << " "
				   << +SEGMENT_LENGTH - 2.0*CORE_RADIUS<< " "
				   << "\n";
	(*outFile) << "-shape -cylinder "
		   << SEGMENT_HEIGHT/2.0 <<" "
		   << ARM_LENGTH/2.0 <<"\n";
	(*outFile) << "-rotation 1 0 0\n";
	(*outFile) << "-mass 0.5 \n";
	(*outFile) << "-colour 0 0.5 0 \n";
	(*outFile) << "-addTouchSensor \n";
	(*outFile) << ") \n";
}

void GA::SendCore(ofstream *outFile) {

	(*outFile) << "\n( Core \n";
	(*outFile) << "-position " << 0.0 << " "
				   << SEGMENT_HEIGHT/2.0  << " "
				   << 0.0 << " "
				   << "\n";
	(*outFile) << "-shape -sphere "
		   << CORE_RADIUS <<"\n";
	(*outFile) << "-mass 1 \n";
	(*outFile) << "-colour 0 1 0 \n";
        (*outFile) << "-addOrientationSensor \n";
	(*outFile) << "-addDistanceSensor TargetObject "<< simParams->targetDistance << " \n";
	(*outFile) << ") \n";
}

void GA::SendEarlyStoppingData(GENOME *currGenome, int currentSlot) {

	char tempFileName[200];
	sprintf(tempFileName,"/tmp/Files%d_%d/temp.dat",simParams->randSeed,currentSlot);

	char fileName[200];
	ofstream *outFile;

	outFile = new ofstream(tempFileName);
		(*outFile) << "1 1\n 0";

	outFile->close();
	delete outFile;
	outFile = NULL;
	sprintf(fileName,"/tmp/Files%d_%d/coreDists.dat",simParams->randSeed,currentSlot);
	simParams->FileRename(tempFileName,fileName);



	outFile = new ofstream(tempFileName);
		(*outFile) << "1 1\n 0";

	outFile->close();
	delete outFile;
	outFile = NULL;
	sprintf(fileName,"/tmp/Files%d_%d/leftShoulderDists.dat",simParams->randSeed,currentSlot);
	simParams->FileRename(tempFileName,fileName);



	outFile = new ofstream(tempFileName);
		(*outFile) << "1 1\n 0";

	outFile->close();
	delete outFile;
	outFile = NULL;
	sprintf(fileName,"/tmp/Files%d_%d/rightShoulderDists.dat",simParams->randSeed,currentSlot);
	simParams->FileRename(tempFileName,fileName);
}

void GA::SendHipGirdles(ofstream *outFile) {

	(*outFile) << "\n( LeftHipGirdle \n";
	(*outFile) << "-position " << -SEGMENT_WIDTH/2.0 << " "
				   << SEGMENT_HEIGHT/2.0  << " "
				   << -SEGMENT_LENGTH + 2.0*CORE_RADIUS<< " "
				   << "\n";
	(*outFile) << "-shape -sphere "
		   << CORE_RADIUS <<"\n";
	(*outFile) << "-mass 1 \n";
	(*outFile) << "-colour 0 1 0 \n";
	(*outFile) << ") \n";

	(*outFile) << "\n( RightHipGirdle \n";
	(*outFile) << "-position " << SEGMENT_WIDTH/2.0 << " "
				   << SEGMENT_HEIGHT/2.0  << " "
				   << -SEGMENT_LENGTH + 2.0*CORE_RADIUS<< " "
				   << "\n";
	(*outFile) << "-shape -sphere "
		   << CORE_RADIUS <<"\n";
	(*outFile) << "-mass 1 \n";
	(*outFile) << "-colour 0 0.5 0 \n";
	(*outFile) << ") \n";
}

void GA::SendHips(ofstream *outFile) {

	(*outFile) << "\n( LeftHip \n";
	(*outFile) << "-position " << -SEGMENT_WIDTH/2.0 << " "
				   << SEGMENT_HEIGHT/2.0  << " "
				   << -SEGMENT_LENGTH + 2.0*CORE_RADIUS<< " "
				   << "\n";
	(*outFile) << "-shape -sphere "
		   << CORE_RADIUS <<"\n";
	(*outFile) << "-mass 1 \n";
	(*outFile) << "-colour 0 1 0 \n";
	(*outFile) << ") \n";

	(*outFile) << "\n( RightHip \n";
	(*outFile) << "-position " << SEGMENT_WIDTH/2.0 << " "
				   << SEGMENT_HEIGHT/2.0  << " "
				   << -SEGMENT_LENGTH + 2.0*CORE_RADIUS<< " "
				   << "\n";
	(*outFile) << "-shape -sphere "
		   << CORE_RADIUS <<"\n";
	(*outFile) << "-mass 1 \n";
	(*outFile) << "-colour 0 0.5 0 \n";
	(*outFile) << ") \n";
}

void GA::SendLeg_Joints(ofstream *outFile) {

	(*outFile) << "\n( BackSegment_LeftHipGirdle\n";
	(*outFile) << "-connect BackSegment LeftHipGirdle\n";
	(*outFile) << "-jointPosition "	<< -SEGMENT_WIDTH/2.0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< -SEGMENT_LENGTH + 2.0*CORE_RADIUS<<" \n";
	(*outFile) << "-jointType Prismatic\n";
	(*outFile) << "-jointNormal -1 0 0\n";
	(*outFile) << "-jointLimits "	<< startRetraction << " "
					<< endExtension << " "
					<< endStanceTransitionWhen << "\n";
	(*outFile) << ")\n";

	(*outFile) << "\n( BackSegment_RightHipGirdle\n";
	(*outFile) << "-connect BackSegment RightHipGirdle\n";
	(*outFile) << "-jointPosition "	<< SEGMENT_WIDTH/2.0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< -SEGMENT_LENGTH + 2.0*CORE_RADIUS<<" \n";
	(*outFile) << "-jointType Prismatic\n";
	(*outFile) << "-jointNormal 1 0 0\n";
	(*outFile) << "-jointLimits "	<< startRetraction << " "
					<< endExtension << " "
					<< endStanceTransitionWhen << "\n";
	(*outFile) << ")\n";


	(*outFile) << "\n( LeftHipGirdle_LeftHip\n";
	(*outFile) << "-connect LeftHipGirdle LeftHip\n";
	(*outFile) << "-jointPosition "	<< -SEGMENT_WIDTH/2.0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< -SEGMENT_LENGTH + 2.0*CORE_RADIUS<<" \n";
	(*outFile) << "-jointType Hinge\n";
	(*outFile) << "-jointNormal 0 0 1\n";
	(*outFile) << "-jointLimits "	<<	stance->Get(2,0)*ratio
					<<" "<<	stance->Get(3,0)*ratio
					<<" "<<	endStance->Get(2,0)*endRatio
					<<" "<<	endStance->Get(3,0)*endRatio
					<<" "<< endStanceTransitionWhen<<" \n";
	(*outFile) << "-addMotor\n";
	(*outFile) << "-motorForce "<<MOTOR_FORCE<<"\n";
	(*outFile) << "-motorSpeed "<<MOTOR_SPEED<<"\n";
        (*outFile) << "-addSensor \n";
	(*outFile) << ")\n";

	(*outFile) << "\n( RightHipGirdle_RightHip\n";
	(*outFile) << "-connect RightHipGirdle RightHip\n";
	(*outFile) << "-jointPosition "	<< SEGMENT_WIDTH/2.0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< -SEGMENT_LENGTH + 2.0*CORE_RADIUS<<" \n";
	(*outFile) << "-jointType Hinge\n";
	(*outFile) << "-jointNormal 0 0 -1\n";
	(*outFile) << "-jointLimits "	<<	stance->Get(2,0)*ratio
					<<" "<<	stance->Get(3,0)*ratio
					<<" "<<	endStance->Get(2,0)*endRatio
					<<" "<<	endStance->Get(3,0)*endRatio
					<<" "<< endStanceTransitionWhen<<" \n";
	(*outFile) << "-addMotor\n";
	(*outFile) << "-motorForce "<<MOTOR_FORCE<<"\n";
	(*outFile) << "-motorSpeed "<<MOTOR_SPEED<<"\n";
        (*outFile) << "-addSensor \n";
	(*outFile) << ")\n";

	(*outFile) << "\n( LeftHip_LeftLeg\n";
	(*outFile) << "-connect LeftHip LeftLeg\n";
	(*outFile) << "-jointPosition "	<< -SEGMENT_WIDTH/2.0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< -SEGMENT_LENGTH + 2.0*CORE_RADIUS<<" \n";
	(*outFile) << "-jointType Hinge\n";
	(*outFile) << "-jointNormal 0 1 0\n";
	(*outFile) << "-jointLimits "	<<	stance->Get(2,1)*ratio
					<<" "<<	stance->Get(3,1)*ratio
					<<" "<<	endStance->Get(2,1)*endRatio
					<<" "<<	endStance->Get(3,1)*endRatio
					<<" "<< endStanceTransitionWhen<<" \n";
	(*outFile) << "-addMotor\n";
	(*outFile) << "-motorForce "<<MOTOR_FORCE<<"\n";
	(*outFile) << "-motorSpeed "<<MOTOR_SPEED<<"\n";
        (*outFile) << "-addSensor \n";
	(*outFile) << ")\n";

	(*outFile) << "\n( RightHip_RightLeg\n";
	(*outFile) << "-connect RightHip RightLeg\n";
	(*outFile) << "-jointPosition "	<< SEGMENT_WIDTH/2.0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< -SEGMENT_LENGTH + 2.0*CORE_RADIUS<<" \n";
	(*outFile) << "-jointType Hinge\n";
	(*outFile) << "-jointNormal 0 -1 0\n";
	(*outFile) << "-jointLimits "	<<	stance->Get(2,1)*ratio
					<<" "<<	stance->Get(3,1)*ratio
					<<" "<<	endStance->Get(2,1)*endRatio
					<<" "<<	endStance->Get(3,1)*endRatio
					<<" "<< endStanceTransitionWhen<<" \n";
	(*outFile) << "-addMotor\n";
	(*outFile) << "-motorForce "<<MOTOR_FORCE<<"\n";
	(*outFile) << "-motorSpeed "<<MOTOR_SPEED<<"\n";
        (*outFile) << "-addSensor \n";
	(*outFile) << ")\n";


	(*outFile) << "\n( LeftLeg_LeftLowerLeg\n";
	(*outFile) << "-connect LeftLeg LeftLowerLeg\n";
	(*outFile) << "-jointPosition "	<< -SEGMENT_WIDTH/2.0 -ARM_LENGTH/2.0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< -SEGMENT_LENGTH + 2.0*CORE_RADIUS<<" \n";
	(*outFile) << "-jointType Prismatic\n";
	(*outFile) << "-jointNormal -1 0 0\n";
	(*outFile) << "-jointLimits "	<< startRetraction << " "
					<< endExtension << " "
					<< endStanceTransitionWhen << "\n";
	(*outFile) << ")\n";

	(*outFile) << "\n( RightLeg_RightLowerLeg\n";
	(*outFile) << "-connect RightLeg RightLowerLeg\n";
	(*outFile) << "-jointPosition "	<< SEGMENT_WIDTH/2.0 +ARM_LENGTH/2.0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< -SEGMENT_LENGTH + 2.0*CORE_RADIUS<<" \n";
	(*outFile) << "-jointType Prismatic\n";
	(*outFile) << "-jointNormal 1 0 0\n";
	(*outFile) << "-jointLimits "	<< startRetraction << " "
					<< endExtension << " "
					<< endStanceTransitionWhen << "\n";
	(*outFile) << ")\n";
}

void GA::SendLegs(ofstream *outFile) {

	(*outFile) << "\n( LeftLeg \n";
	(*outFile) << "-position " << -SEGMENT_WIDTH/2.0-ARM_LENGTH/4.0 << " "
				   << SEGMENT_HEIGHT/2.0 << " "
				   << -SEGMENT_LENGTH + 2.0*CORE_RADIUS<< " "
				   << "\n";
	(*outFile) << "-shape -cylinder "
		   << SEGMENT_HEIGHT/2.0 <<" "
		   << ARM_LENGTH/2.0 <<"\n";
	(*outFile) << "-rotation 1 0 0\n";
	(*outFile) << "-mass 0.5 \n";
	(*outFile) << "-colour 0 1 0 \n";
	(*outFile) << ") \n";

	(*outFile) << "\n( RightLeg \n";
	(*outFile) << "-position " << +SEGMENT_WIDTH/2.0+ARM_LENGTH/4.0 << " "
				   << SEGMENT_HEIGHT/2.0 << " "
				   << -SEGMENT_LENGTH + 2.0*CORE_RADIUS<< " "
				   << "\n";
	(*outFile) << "-shape -cylinder "
		   << SEGMENT_HEIGHT/2.0 <<" "
		   << ARM_LENGTH/2.0 <<"\n";
	(*outFile) << "-rotation 1 0 0\n";
	(*outFile) << "-mass 0.5 \n";
	(*outFile) << "-colour 0 0.5 0 \n";
	(*outFile) << ") \n";


	(*outFile) << "\n( LeftLowerLeg \n";
	(*outFile) << "-position " << -SEGMENT_WIDTH/2.0-3.0*ARM_LENGTH/4.0 << " "
				   << SEGMENT_HEIGHT/2.0 << " "
				   << -SEGMENT_LENGTH + 2.0*CORE_RADIUS<< " "
				   << "\n";
	(*outFile) << "-shape -cylinder "
		   << SEGMENT_HEIGHT/2.0 <<" "
		   << ARM_LENGTH/2.0 <<"\n";
	(*outFile) << "-rotation 1 0 0\n";
	(*outFile) << "-mass 0.5 \n";
	(*outFile) << "-colour 0 1 0 \n";
	(*outFile) << "-addTouchSensor \n";
	(*outFile) << ") \n";

	(*outFile) << "\n( RightLowerLeg \n";
	(*outFile) << "-position " << +SEGMENT_WIDTH/2.0+3.0*ARM_LENGTH/4.0 << " "
				   << SEGMENT_HEIGHT/2.0 << " "
				   << -SEGMENT_LENGTH + 2.0*CORE_RADIUS<< " "
				   << "\n";
	(*outFile) << "-shape -cylinder "
		   << SEGMENT_HEIGHT/2.0 <<" "
		   << ARM_LENGTH/2.0 <<"\n";
	(*outFile) << "-rotation 1 0 0\n";
	(*outFile) << "-mass 0.5 \n";
	(*outFile) << "-colour 0 0.5 0 \n";
	(*outFile) << "-addTouchSensor \n";
	(*outFile) << ") \n";
}

void GA::SendLowerArm(ofstream *outFile) {

	(*outFile) << "\n( LowerArm \n";
	(*outFile) << "-position " << "0 "
                                   << ARM_HEIGHT  << " "
                                   << UPPERARM_LENGTH + LOWERARM_LENGTH/2.0 << " "
				   << "\n";
	(*outFile) << "-shape -cylinder "	<<ARM_RADIUS<<" "
						<<LOWERARM_LENGTH<<"\n";
	(*outFile) << "-rotation 0 0 1\n";
	(*outFile) << "-mass 1 \n";
	(*outFile) << "-colour 0.5 0.5 0.5 \n";
	(*outFile) << ") \n";
}

void GA::SendLowerArm_Joint(ofstream *outFile) {

	(*outFile) << "\n( LowerArm_UpperArm\n";
	(*outFile) << "-connect LowerArm UpperArm\n";
	(*outFile) << "-jointPosition "	<< 0 << " "
				 	<< ARM_HEIGHT << " "
				     	<< UPPERARM_LENGTH <<" \n";
	(*outFile) << "-jointType Hinge\n";
	(*outFile) << "-jointNormal 0 1 0\n";
	(*outFile) << "-jointLimits "<<-LOWERARM_JOINT_RANGE<<" "<<LOWERARM_JOINT_RANGE<<" \n";
	(*outFile) << "-addSensor\n";
	(*outFile) << "-addMotor\n";
	(*outFile) << "-motorForce "<<ARM_MOTOR_FORCE<<"\n";
	(*outFile) << "-motorSpeed "<<MOTOR_SPEED<<"\n";
	(*outFile) << ")\n";
}

void GA::SendPostscript(ofstream *outFile) {

	(*outFile) << "\n( \n";
	(*outFile) << "-evaluationPeriod "<< simParams->evaluationLength <<" \n";
	(*outFile) << "-testForExplosions \n";
	(*outFile) << ") \n";
}

void GA::SendSegments(ofstream *outFile) {

	(*outFile) << "\n( FrontSegment \n";
	(*outFile) << "-position " << 0.0 << " "
				   << SEGMENT_HEIGHT/2.0  << " "
				   << SEGMENT_LENGTH/2.0 << " "
				   << "\n";
	(*outFile) << "-shape -rectangle "
		   << SEGMENT_LENGTH <<" "
		   << SEGMENT_WIDTH <<" "
		   << SEGMENT_HEIGHT <<"\n";
	(*outFile) << "-mass 1 \n";
	(*outFile) << "-colour 1 0 0 \n";
	(*outFile) << ") \n";

	(*outFile) << "\n( BackSegment \n";
	(*outFile) << "-position " << 0.0 << " "
				   << SEGMENT_HEIGHT/2.0  << " "
				   << -SEGMENT_LENGTH/2.0 << " "
				   << "\n";
	(*outFile) << "-shape -rectangle "
		   << SEGMENT_LENGTH <<" "
		   << SEGMENT_WIDTH <<" "
		   << SEGMENT_HEIGHT <<"\n";
	(*outFile) << "-mass 1 \n";
	(*outFile) << "-colour 0 0 1 \n";
	(*outFile) << ") \n";
}

void GA::SendShoulder(ofstream *outFile) {

	(*outFile) << "\n( Shoulder \n";
	(*outFile) << "-position " << "0 "
                                   << ARM_HEIGHT  << " "
                                   << 0 << " "
				   << "\n";
	(*outFile) << "-shape -sphere "	<<SHOULDER_RADIUS<<"\n";
	(*outFile) << "-mass 1 \n";
	(*outFile) << "-colour 0 0 1 \n";
	(*outFile) << ") \n";
}

void GA::SendShoulderGirdles(ofstream *outFile) {

	(*outFile) << "\n( LeftShoulderGirdle \n";
	(*outFile) << "-position " << -SEGMENT_WIDTH/2.0 << " "
				   << SEGMENT_HEIGHT/2.0  << " "
				   << +SEGMENT_LENGTH - 2.0*CORE_RADIUS<< " "
				   << "\n";
	(*outFile) << "-shape -sphere "
		   << CORE_RADIUS <<"\n";
	(*outFile) << "-mass 1 \n";
	(*outFile) << "-colour 0 1 0 \n";
	(*outFile) << ") \n";

	(*outFile) << "\n( RightShoulderGirdle \n";
	(*outFile) << "-position " << SEGMENT_WIDTH/2.0 << " "
				   << SEGMENT_HEIGHT/2.0  << " "
				   << +SEGMENT_LENGTH - 2.0*CORE_RADIUS<< " "
				   << "\n";
	(*outFile) << "-shape -sphere "
		   << CORE_RADIUS <<"\n";
	(*outFile) << "-mass 1 \n";
	(*outFile) << "-colour 0 0.5 0 \n";
	(*outFile) << ") \n";
}

void GA::SendShoulders(ofstream *outFile) {

	(*outFile) << "\n( LeftShoulder \n";
	(*outFile) << "-position " << -SEGMENT_WIDTH/2.0 << " "
				   << SEGMENT_HEIGHT/2.0  << " "
				   << SEGMENT_LENGTH - 2.0*CORE_RADIUS << " "
				   << "\n";
	(*outFile) << "-shape -sphere "
		   << CORE_RADIUS <<"\n";
	(*outFile) << "-mass 1 \n";
	(*outFile) << "-colour 0 1 0 \n";
	(*outFile) << "-addDistanceSensor TargetObject "<< simParams->targetDistance << " \n";
	(*outFile) << ") \n";

	(*outFile) << "\n( RightShoulder \n";
	(*outFile) << "-position " << SEGMENT_WIDTH/2.0 << " "
				   << SEGMENT_HEIGHT/2.0  << " "
				   << SEGMENT_LENGTH - 2.0*CORE_RADIUS << " "
				   << "\n";
	(*outFile) << "-shape -sphere "
		   << CORE_RADIUS <<"\n";
	(*outFile) << "-mass 1 \n";
	(*outFile) << "-colour 0 0.5 0 \n";
	(*outFile) << "-addDistanceSensor TargetObject "<< simParams->targetDistance << " \n";
	(*outFile) << ") \n";
}

void GA::SendSpine_Joints(ofstream *outFile) {

	(*outFile) << "\n( Core_FrontSegment\n";
	(*outFile) << "-connect Core FrontSegment\n";
	(*outFile) << "-jointPosition "	<< 0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< 0 <<" \n";
	(*outFile) << "-jointType Hinge\n";
	(*outFile) << "-jointNormal 0 0 1\n";
	(*outFile) << "-jointLimits "	<<-SPINE_VERTICAL_JOINT_RANGE<<" "<<SPINE_VERTICAL_JOINT_RANGE<<" "
					<<-SPINE_VERTICAL_JOINT_RANGE<<" "<<SPINE_VERTICAL_JOINT_RANGE<<" 1.0 \n";
	(*outFile) << "-addMotor\n";
	(*outFile) << "-motorForce "<<MOTOR_FORCE<<"\n";
	(*outFile) << "-motorSpeed "<<MOTOR_SPEED<<"\n";
        (*outFile) << "-addSensor \n";
	(*outFile) << ")\n";

	(*outFile) << "\n( Core_BackSegment\n";
	(*outFile) << "-connect Core BackSegment\n";
	(*outFile) << "-jointPosition "	<< 0 << " "
				 	<< SEGMENT_HEIGHT/2.0 << " "
				     	<< 0 <<" \n";
	(*outFile) << "-jointType Hinge\n";
	(*outFile) << "-jointNormal 0 1 0\n";
	(*outFile) << "-jointLimits "	<<-SPINE_HORIZONTAL_JOINT_RANGE<<" "<<SPINE_HORIZONTAL_JOINT_RANGE<<" "
					<<-SPINE_HORIZONTAL_JOINT_RANGE<<" "<<SPINE_HORIZONTAL_JOINT_RANGE<<" 1.0 \n";
	(*outFile) << "-addMotor\n";
	(*outFile) << "-motorForce "<<MOTOR_FORCE<<"\n";
	(*outFile) << "-motorSpeed "<<MOTOR_SPEED<<"\n";
        (*outFile) << "-addSensor \n";
	(*outFile) << ")\n";
}

void GA::SendTargetObject(ofstream *outFile, int objectIndex) {

	(*outFile) << "\n( TargetObject \n";

//	double theta   = -3.14159/2.0;
	double theta = simParams->Scale(double(objectIndex),
					0.0,double(simParams->totalObjects-1),
					-simParams->turningAmount,simParams->turningAmount);
	double xOffset = simParams->targetDistance * sin(theta);
	double yOffset = simParams->targetDistance * cos(theta);

	(*outFile) << "-position " << xOffset << " "
				   << TARGET_LWH/2.0  << " "
				   << yOffset << " "
				   << "\n";

	(*outFile) << "-shape -rectangle "
		   << TARGET_LWH <<" "
		   << TARGET_LWH <<" "
		   << TARGET_LWH <<"\n";

	(*outFile) << "-mass 1 \n";
	(*outFile) << "-colour 0.2 0.2 0.2 \n";
	(*outFile) << ") \n";
}

void GA::SendUpperArm(ofstream *outFile) {

	(*outFile) << "\n( UpperArm \n";
	(*outFile) << "-position " << "0 "
                                   << ARM_HEIGHT  << " "
                                   << UPPERARM_LENGTH/2.0 << " "
				   << "\n";
	(*outFile) << "-shape -cylinder "	<<ARM_RADIUS<<" "
						<<UPPERARM_LENGTH<<"\n";
	(*outFile) << "-rotation 0 0 1\n";
	(*outFile) << "-mass 1 \n";
	(*outFile) << "-colour 0.5 0.5 0.5 \n";
	(*outFile) << ") \n";
}

void GA::SendUpperArm_Joint(ofstream *outFile) {

	(*outFile) << "\n( UpperArm_Shoulder\n";
	(*outFile) << "-connect UpperArm Shoulder\n";
	(*outFile) << "-jointPosition "	<< 0 << " "
				 	<< ARM_HEIGHT << " "
				     	<< 0 <<" \n";
	(*outFile) << "-jointType Hinge\n";
	(*outFile) << "-jointNormal 1 0 0\n";
	(*outFile) << "-jointLimits "<<-UPPERARM_JOINT_RANGE<<" "<<UPPERARM_JOINT_RANGE<<" \n";
	(*outFile) << "-addSensor\n";
	(*outFile) << "-addMotor\n";
	(*outFile) << "-motorForce "<<ARM_MOTOR_FORCE<<"\n";
	(*outFile) << "-motorSpeed "<<MOTOR_SPEED<<"\n";
	(*outFile) << ")\n";
}

void GA::SeverChildConnections(GENOME *dyingParent) {

	for (int i=0;i<simParams->popSize;i++) {

		if ( genomes[i] )
			if ( genomes[i]!=dyingParent )
				if ( genomes[i]->myParent == dyingParent )
					genomes[i]->myParent = NULL;
	}
}

int GA::SlotFinished(int currentSlot) {

	int finished = true;

	char fileName[200];

	sprintf(fileName,"/tmp/Files%d_%d/myCoreDists.dat",simParams->randSeed,currentSlot);
	if ( finished && (!simParams->FileExists(fileName)) )
		finished = false;

	sprintf(fileName,"/tmp/Files%d_%d/myLeftShoulderDists.dat",simParams->randSeed,currentSlot);
	if ( finished && (!simParams->FileExists(fileName)) )
		finished = false;

	sprintf(fileName,"/tmp/Files%d_%d/myRightShoulderDists.dat",simParams->randSeed,currentSlot);
	if ( finished && (!simParams->FileExists(fileName)) )
		finished = false;

	sprintf(fileName,"/tmp/Files%d_%d/Sensors.dat",simParams->randSeed,currentSlot);
	if ( finished && (!simParams->FileExists(fileName)) )
		finished = false;

	return( finished );
}

void GA::SortGenomesOnLayer(int layerIndex) {

        int firstIndexInLayer = simParams->numIndsPerLayer * (layerIndex);
        int lastIndexInLayer  = simParams->numIndsPerLayer * (layerIndex+1) - 1;

        for (int i=firstIndexInLayer;i<=lastIndexInLayer-1;i++) {

		for (int j=i+1;j<=lastIndexInLayer;j++) {

			if ( genomes[j]->fitness > genomes[i]->fitness ) {
				
				GENOME *temp = genomes[i];
				genomes[i] = genomes[j];
				genomes[j] = temp;
				temp = NULL;
			}
		}
	}
}

void GA::SortGenomesOnLayers(void) {

	for (int i=0;i<simParams->numLayers;i++)

		SortGenomesOnLayer(i);
}

void GA::Stance_Create(void) {

	stance_Prone = new MATRIX(4,4,0.0);					
	stance_Prone->Set(0,0,0		-LEG_VERTICAL_JOINT_RANGE);	// Spine to shoulder	
	stance_Prone->Set(0,1,0		-LEG_HORIZONTAL_JOINT_RANGE); 	// shoulder to upper leg	
	stance_Prone->Set(0,2,0		-LEG_HORIZONTAL_JOINT_RANGE); 	// upper leg to middle leg
	stance_Prone->Set(0,3,0		-LEG_HORIZONTAL_JOINT_RANGE);	// middle leg to lower leg
	
	stance_Prone->Set(1,0,0		+LEG_VERTICAL_JOINT_RANGE); 	
	stance_Prone->Set(1,1,0		+LEG_HORIZONTAL_JOINT_RANGE);
 	stance_Prone->Set(1,2,0		+LEG_HORIZONTAL_JOINT_RANGE);
 	stance_Prone->Set(1,3,0		+LEG_HORIZONTAL_JOINT_RANGE);

	stance_Prone->Set(2,0,0		-LEG_VERTICAL_JOINT_RANGE);
 	stance_Prone->Set(2,1,0		-LEG_HORIZONTAL_JOINT_RANGE);
 	stance_Prone->Set(2,2,0		-LEG_HORIZONTAL_JOINT_RANGE);
 	stance_Prone->Set(2,3,0		-LEG_HORIZONTAL_JOINT_RANGE);

	stance_Prone->Set(3,0,0		+LEG_VERTICAL_JOINT_RANGE);
 	stance_Prone->Set(3,1,0		+LEG_HORIZONTAL_JOINT_RANGE);
 	stance_Prone->Set(3,2,0		+LEG_HORIZONTAL_JOINT_RANGE);
 	stance_Prone->Set(3,3,0		+LEG_HORIZONTAL_JOINT_RANGE);

	stance_Upright = new MATRIX(stance_Prone);					
	stance_Upright->Set(0,0,stance_Upright->Get(0,0)-90);	
	stance_Upright->Set(1,0,stance_Upright->Get(1,0)-90); 	
	stance_Upright->Set(2,0,stance_Upright->Get(2,0)-90);
	stance_Upright->Set(3,0,stance_Upright->Get(3,0)-90);
	// Bend legs to vertical

	double theta = 0;
	stance_Four_Feet = new MATRIX(stance_Upright);
	stance_Four_Feet->Set(0,1,stance_Four_Feet->Get(0,1)+theta);	
	stance_Four_Feet->Set(1,1,stance_Four_Feet->Get(1,1)+theta);
	stance_Four_Feet->Set(0,2,stance_Four_Feet->Get(0,2)-2.0*theta);	
	stance_Four_Feet->Set(1,2,stance_Four_Feet->Get(1,2)-2.0*theta);	
	stance_Four_Feet->Set(0,3,stance_Four_Feet->Get(0,3)+(90+theta));	
	stance_Four_Feet->Set(1,3,stance_Four_Feet->Get(1,3)+(90+theta));

	stance_Four_Feet->Set(2,1,stance_Four_Feet->Get(2,1)+theta);	
	stance_Four_Feet->Set(3,1,stance_Four_Feet->Get(3,1)+theta);
	stance_Four_Feet->Set(2,2,stance_Four_Feet->Get(2,2)-2.0*theta);	
	stance_Four_Feet->Set(3,2,stance_Four_Feet->Get(3,2)-2.0*theta);	
	stance_Four_Feet->Set(2,3,stance_Four_Feet->Get(2,3)+(90+theta));	
	stance_Four_Feet->Set(3,3,stance_Four_Feet->Get(3,3)+(90+theta));
	// Create four feet

	stance_Squat_Upright = new MATRIX(stance_Upright);
	stance_Squat_Upright->Set(0,1,stance_Squat_Upright->Get(0,1)-45);	
	stance_Squat_Upright->Set(1,1,stance_Squat_Upright->Get(1,1)-45);
	stance_Squat_Upright->Set(0,2,stance_Squat_Upright->Get(0,2)+90);	
	stance_Squat_Upright->Set(1,2,stance_Squat_Upright->Get(1,2)+90);	
	stance_Squat_Upright->Set(0,3,stance_Squat_Upright->Get(0,3)+45);	
	stance_Squat_Upright->Set(1,3,stance_Squat_Upright->Get(1,3)+45);	

	stance_Squat_Upright->Set(2,1,stance_Squat_Upright->Get(2,1)+80);	
	stance_Squat_Upright->Set(3,1,stance_Squat_Upright->Get(3,1)+80);
	stance_Squat_Upright->Set(2,2,stance_Squat_Upright->Get(2,2)-160);	
	stance_Squat_Upright->Set(3,2,stance_Squat_Upright->Get(3,2)-160);
	stance_Squat_Upright->Set(2,3,stance_Squat_Upright->Get(2,3)+160);	
	stance_Squat_Upright->Set(3,3,stance_Squat_Upright->Get(3,3)+160);
	// Accordian legs inward

	stance_Squat = new MATRIX(stance_Squat_Upright);
	stance_Squat->Set(0,0,stance_Squat->Get(0,0)+90);	
	stance_Squat->Set(1,0,stance_Squat->Get(1,0)+90); 	
	stance_Squat->Set(2,0,stance_Squat->Get(2,0)+90);
	stance_Squat->Set(3,0,stance_Squat->Get(3,0)+90);
	// Bend legs back to horizontal

	stance = new MATRIX(stance_Prone);

	endStance = new MATRIX(stance_Four_Feet);

	endStanceTransitionWhen = 1.0;

	startRetraction = -ARM_LENGTH/2.0;
	endExtension  = -ARM_LENGTH/2.0;

	ratio = 1.0;
	endRatio = 1.0;

}

void GA::Stance_Destroy(void) {

	delete endStance;
	endStance = NULL;

	delete stance;
	stance = NULL;

	delete stance_Four_Feet;
	stance_Four_Feet = NULL;

	delete stance_Prone;
	stance_Prone = NULL;

	delete stance_Upright;
	stance_Upright = NULL;

	delete stance_Squat;
	stance_Squat = NULL;

	delete stance_Squat_Upright;
	stance_Squat_Upright = NULL;
}

void GA::Stance_SetForBetweenScaffolding(int objectIndex) {

        ratio = double(objectIndex+1) / double(simParams->totalObjects);
        endRatio = ratio;

        delete stance;
        stance = new MATRIX(stance_Four_Feet);
        stance->Mult(ratio);
        MATRIX *temp = new MATRIX(stance_Prone);
        temp->Mult(1.0-ratio);
        stance->Add(temp);
        delete temp;
        temp = NULL;

        delete endStance;
        endStance = new MATRIX(stance_Four_Feet);
        endStance->Mult(endRatio);
        temp = new MATRIX(stance_Prone);
        temp->Mult(1.0-endRatio);
        endStance->Add(temp);
        delete temp;
        temp = NULL;

	startRetraction = -(ARM_LENGTH/2.0) +(ARM_LENGTH/2.0)*ratio;
	endExtension  = -(ARM_LENGTH/2.0) +(ARM_LENGTH/2.0)*endRatio;
}

void GA::Stance_SetForNoScaffolding(void) {

	ratio = 1.0;
	endRatio = 1.0;

        delete stance;
        stance = new MATRIX(stance_Four_Feet);

        delete endStance;
        endStance = new MATRIX(stance_Four_Feet);

	startRetraction = -(ARM_LENGTH/2.0) +(ARM_LENGTH/2.0)*ratio;
	endExtension  = -(ARM_LENGTH/2.0) +(ARM_LENGTH/2.0)*endRatio;
}


void GA::Stance_SetForWithinScaffolding(int objectIndex) {

        double halfTheObjects = double(simParams->totalObjects)/2.0;

        if ( objectIndex<halfTheObjects ) {

		double increaseIncrement = 0.25 + objectIndex*(1.0-0.25)/(halfTheObjects-1.0);

                ratio = 0.25;

		endRatio = increaseIncrement;
        }
        else {

		double objectIndexAfterHalfwayPoint = objectIndex - halfTheObjects;

		double increaseIncrement = 0.25 + (objectIndexAfterHalfwayPoint+1)*(1.0-0.25)/(halfTheObjects);

		ratio = increaseIncrement;
		
                endRatio = 1.0;
        }

        delete stance;
        stance = new MATRIX(stance_Four_Feet);
        stance->Mult(ratio);
        MATRIX *temp = new MATRIX(stance_Prone);
        temp->Mult(1.0-ratio);
        stance->Add(temp);
        delete temp;
        temp = NULL;

        delete endStance;
        endStance = new MATRIX(stance_Four_Feet);
        endStance->Mult(endRatio);
        temp = new MATRIX(stance_Prone);
        temp->Mult(1.0-endRatio);
        endStance->Add(temp);
        delete temp;
        temp = NULL;

	startRetraction = -(ARM_LENGTH/2.0) +(ARM_LENGTH/2.0)*ratio;
	endExtension  = -(ARM_LENGTH/2.0) +(ARM_LENGTH/2.0)*endRatio;
}

void GA::Stance_SetForWithinScaffolding2(int objectIndex) {

        endRatio = 1.0;

        for (int i=0;i<simParams->totalObjects;i++)

		if ( i==objectIndex )
                	ratio = double(i)/(double(simParams->totalObjects)-1.0);

        delete stance;
        stance = new MATRIX(stance_Four_Feet);
        stance->Mult(ratio);
        MATRIX *temp = new MATRIX(stance_Prone);
        temp->Mult(1.0-ratio);
        stance->Add(temp);
        delete temp;
        temp = NULL;

        delete endStance;
        endStance = new MATRIX(stance_Four_Feet);
        endStance->Mult(endRatio);
        temp = new MATRIX(stance_Prone);
        temp->Mult(1.0-endRatio);
        endStance->Add(temp);
        delete temp;
        temp = NULL;

	startRetraction = -(ARM_LENGTH/2.0) +(ARM_LENGTH/2.0)*ratio;
	endExtension  = -(ARM_LENGTH/2.0) +(ARM_LENGTH/2.0)*endRatio;
}

void GA::Stance_SetForWithinScaffolding3(int objectIndex) {

        ratio = 0.0;
        endRatio = 1.0;

	for (int i=0;i<simParams->totalObjects;i++)

		if ( i==objectIndex )
			endStanceTransitionWhen = double(simParams->totalObjects-(i+1))/double(simParams->totalObjects-1);

        delete stance;
        stance = new MATRIX(stance_Four_Feet);
        stance->Mult(ratio);
        MATRIX *temp = new MATRIX(stance_Prone);
        temp->Mult(1.0-ratio);
        stance->Add(temp);
        delete temp;
        temp = NULL;

        delete endStance;
        endStance = new MATRIX(stance_Four_Feet);
        endStance->Mult(endRatio);
        temp = new MATRIX(stance_Prone);
        temp->Mult(1.0-endRatio);
        endStance->Add(temp);
        delete temp;
        temp = NULL;

	startRetraction = -(ARM_LENGTH/2.0) +(ARM_LENGTH/2.0)*ratio;
	endExtension  = -(ARM_LENGTH/2.0) +(ARM_LENGTH/2.0)*endRatio;
}

void GA::Stance_SetForWithinScaffolding4(int objectIndex) {

        double halfTheObjects = double(simParams->totalObjects)/2.0;

        if ( (objectIndex%2)==0 ) {

                ratio    = ((double(objectIndex)/2.0)+0.0) / halfTheObjects;
                endRatio = ((double(objectIndex)/2.0)+1.0) / halfTheObjects;

        }
        else {
                ratio    = ((double(objectIndex-1)/2.0)+1.0) / halfTheObjects;
                endRatio = ((double(objectIndex-1)/2.0)+1.0) / halfTheObjects;
        }

        delete stance;
        stance = new MATRIX(stance_Four_Feet);
        stance->Mult(ratio);
        MATRIX *temp = new MATRIX(stance_Prone);
        temp->Mult(1.0-ratio);
        stance->Add(temp);
        delete temp;
        temp = NULL;

        delete endStance;
        endStance = new MATRIX(stance_Four_Feet);
        endStance->Mult(endRatio);
        temp = new MATRIX(stance_Prone);
        temp->Mult(1.0-endRatio);
        endStance->Add(temp);
        delete temp;
        temp = NULL;

	startRetraction = -(ARM_LENGTH/2.0) +(ARM_LENGTH/2.0)*ratio;
	endExtension  = -(ARM_LENGTH/2.0) +(ARM_LENGTH/2.0)*endRatio;
}

int  GA::SuccessfulGenomePresent(void) {

	int i=0;
	int found = false;

	while ( (i<simParams->popSize) && (!found) ) {

		if ( genomes[i] && genomes[i]->Successful() )
			found = true;
		else
			i++;
	}

	return( found );
}

int GA::SuperiorOnLayer(int genomeIndex, int layerIndex) {

	int superior = false;
	int i = 0;

	while ( 	(!superior) && 
			(i<simParams->popSize) ) {

		if ( 	genomes[i] &&
			(genomes[genomeIndex]!=genomes[i]) &&
			(AgeLayer(i)==layerIndex) &&
			genomes[genomeIndex]->SuperiorTo(genomes[i]) )

			superior = true;
		else
			i++; 

	}

	return( superior );
}

void GA::SwapParentPointers(GENOME *firstParent, GENOME *secondParent) {

	for (int i=0;i<simParams->popSize;i++)
		genomes[i]->SwapParentPointer(firstParent,secondParent);
}

int  GA::TooMuchExperienceForLayer(int genomeIndex, int layerIndex) {

	return( genomes[genomeIndex]->objectsToEvaluate > ObjectsForLayer(layerIndex) );
}

int    GA::TooOld(int genomeIndex) {

	if ( AgeLayer(genomeIndex)==(simParams->numLayers-1) )
		return false;
	
	int tooOld = genomes[genomeIndex]->Age(int(numGenomesEvaluated)) > 
                	MaxAgeForLayer(AgeLayer(genomeIndex));

	return( tooOld );
}

void GA::TryMoveUp(int i, int slotIndex) {

	if ( Empty(i) )
		return;

	if ( AgeLayer(i)==(simParams->numLayers-1) )
		return;

	ReplaceVictimAbove(i,slotIndex);
}

void GA::TurnOffReinitialization(void) {

	reinitializationMode = false;
}

void GA::TurnOnReinitialization(void) {

	reinitializationMode = true;
	reinitializationIndividual = 0;
}

int  GA::VictimFailed(int aggressor, int victimIndex, int slotIndex) {

	if ( Empty(victimIndex) )
		return( true );

	if ( TooOld(victimIndex) )
		return( true );

	if ( TooMuchExperienceForLayer(aggressor,AgeLayer(victimIndex)) )
		return( true );

	if ( genomes[aggressor]->Dominates(genomes[victimIndex]) )

		return( HandleAggressorDominatingVictim(aggressor,victimIndex,slotIndex) );


	if ( genomes[victimIndex]->Dominates(genomes[aggressor]) )

		return( false );

	return( HandleNeitherDominatingOneAnother(aggressor, victimIndex) );
}

void GA::WaitForSensorReadings(void) {

/*
	simParams->WaitForFile(simParams->sensorFileName);

	char fileName[200];

	sprintf(fileName,"/tmp/Files%d/hiddenNeuronReadings.dat",simParams->randSeed);
	simParams->WaitForFile(fileName);

	sprintf(fileName,"/tmp/Files%d/myAPs.dat",simParams->randSeed);
	simParams->WaitForFile(fileName);
*/
}

void GA::WriteBest(ofstream *outFile, int layerIndex) {

	int bestIndex = -1;
	double maxFit = -1000.0;
	for (int i=0;i<simParams->popSize;i++) {
		if ( 	genomes[i] && 
			(AgeLayer(i)==layerIndex) &&
			((!genomes[i]->dominated)||(NumNonDominated(layerIndex)==0)) ) {
			if ( genomes[i]->fitness > maxFit ) {
				maxFit = genomes[i]->fitness;
				bestIndex = i;
			}
		}
	}

	if ( bestIndex==-1 )
		currBestFit = 	0;
	else
		currBestFit = 	genomes[bestIndex]->coreDist * 
				genomes[bestIndex]->leftShoulderDist * 
				genomes[bestIndex]->rightShoulderDist;

	(*outFile) << currBestFit << " ";
}

void GA::WriteMean(ofstream *outFile, int layerIndex) {

	int bestIndex;
	double meanFit = 0.0;
	double nums = 0.0;

	for (int i=0;i<simParams->popSize;i++) {

		if ( 	genomes[i] && 
			(AgeLayer(i)==layerIndex) &&
			((!genomes[i]->dominated)||(NumNonDominated(layerIndex)==0)) ) {
			meanFit = meanFit + (
					genomes[i]->coreDist*
					genomes[i]->leftShoulderDist*
					genomes[i]->rightShoulderDist);
			nums++;
		}
	}
	
	if ( nums>0.0 )
		meanFit = meanFit / nums;

	(*outFile) << meanFit << " ";
}

void GA::WriteAll(void) {

	SaveBest();

	ofstream *outFile = new ofstream(dataFileName,ios::app);

//	(*outFile) << MostObjectsBeingTackled() << " ";

	(*outFile) << simParams->stanceIndex << " ";

	(*outFile) << numGenomesEvaluated 			<< " ";

	double cpuMinsSinceStart = 	(clock()/(60.0*CLOCKS_PER_SEC)) - 
					(simParams->startingClock/(60.0*CLOCKS_PER_SEC));
	(*outFile) << cpuMinsSinceStart 			<< " ";

	for (int i=simParams->numLayers-1;i>=0;i=i-1)
		(*outFile) << ObjectsForLayer(i) << " ";

	for (int i=simParams->numLayers-1;i>=0;i=i-1)
		WriteBest(outFile,i);

	for (int i=simParams->numLayers-1;i>=0;i=i-1)
		WriteMean(outFile,i);

	for (int i=simParams->numLayers-1;i>=0;i=i-1)
		(*outFile) <<		NumNonDominated(i)	<< " ";

	(*outFile) << "\n";

	outFile->close();
	delete outFile;

	outputDataWritten	 = true;
	meanObjectsEvaluated 	= 0.0;
	meanEvalTimes 		= 0.0;
	meanCoreDists		= 0.0;
	meanLeftShoulderDists			= 0.0;
	meanRightShoulderDists		= 0.0;
	bestCoreDist		= 0.0;
	bestLeftShoulderDist			= 0.0;
	bestRightShoulderDist		= 0.0;
	nums 			= 0.0;
}

int    GA::YoungEnoughParent(int genomeIndex) {

	return( genomes[genomeIndex]->Age(int(numGenomesEvaluated)) < 
                MaxAgeForLayer(AgeLayer(genomeIndex)) );
}

#endif
