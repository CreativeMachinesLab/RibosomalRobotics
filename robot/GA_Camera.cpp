// GA_Camera.cpp : Defines the entry point for the console application.
//

#include "ga.h"
#include "simParams.h"

//#ifdef _DEBUG
//#define new DEBUG_NEW
//#undef THIS_FILE
//static char THIS_FILE[] = __FILE__;
//#endif

SIM_PARAMS *simParams;

extern char FILES_LOCATION[100];
extern char DATA_DIRECTORY[100];

int main(int argc, char* argv[], char* envp[])
{
	simParams = new SIM_PARAMS(argc,argv);

	GA *ga = new GA(simParams->loadState);
	
	ga->Evolve();

	if ( simParams->TimeElapsed() )
		ga->SaveState();

	delete ga;
	ga  = NULL;

	delete simParams;
	simParams = NULL;

	return( 0 );
}


