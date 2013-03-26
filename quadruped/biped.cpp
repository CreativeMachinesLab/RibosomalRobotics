//QUADRUPED

/*************************************************************************
*                                                                       *
* Open Dynamics Engine, Copyright (C) 2001,2002 Russell L. Smith.       *
* All rights reserved.  Email: russ@q12.org   Web: www.q12.org          *
*                                                                       *R
* This library is free software; you can redistribute it and/or         *
* modify it under the terms of EITHER:                                  *
*   (1) The GNU Lesser General Public License as published by the Free  *
*       Software Foundation; either version 2.1 of the License, or (at  *
*       your option) any later version. The text of the GNU Lesser      *
*       General Public License is included with this library in the     *
*       file LICENSE.TXT.                                               *
*   (2) The BSD-style license that is included with this library in     *
*       the file LICENSE-BSD.TXT.                                       *
*                                                                       *
* This library is distributed in the hope that it will be useful,       *
* but WITHOUT ANY WARRANTY; without even the implied warranty of        *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the files    *
* LICENSE.TXT and LICENSE-BSD.TXT for more details.                     *
*                                                                       *
*************************************************************************/
//evolve output_ fitness bipedstartgenes_leo substrate2i2df_quad 0 4

//http://comments.gmane.org/gmane.comp.lib.ode/9510

#define GRAPHICS

#include <ode/odeconfig.h>
#include "ode/ode.h"
#include <assert.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h> 
#include <time.h>
#include "Controller.h"
#include "Quadruped.h"
#include <sys/types.h>
//#include <dirent.h>
#include <errno.h>

#ifdef GRAPHICS
	#include <tchar.h> 
#endif

#include <float.h>
#include <string.h>
#include <cstring>

//#include <unistd.h>


#include "neat.h"
#include "organism.h"
#include "noveltyset.h"
#include "datarec.h"

#include <sys/stat.h>
//#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef GRAPHICS
	#include <windows.h>
#endif

using namespace std;

//#include "ga.h"
//#include "ConfigFile.h"
#include <iostream>
//0.25
static bool WATER=false;

static int NUMBER_SEGMENTS =2; //0.2 - 0.8
static bool RESET=false;
static int SHAPE_MODE=2; //2 = robots
static bool STEP_MODE=false;

//these variables are used if the CPPN determines the number of segments
static bool CPPN_DET_NUM_PARTS=false;
static int SEGMENTS_BETWEEN_MOTORS; 

//0.0, 0.7, 1.4
static double LEG_ANGEL =0.0; //0.3,  0.6 1.0

static dReal MAXTORQUE_HIPMINOR= 5.0;  // 5.0
static dReal MAXTORQUE_HIPMAJOR= 5; //5.0

static bool MOTORS_ENABLED = true;

static int const CPPN_QUERIES = 10;

static float const LINK_PROB = 0.7;
static int const CPPN_INPUTS = 2; //was 3
static int const CPPN_HIDDEN = 2;
static int CPPN_OUTPUTS = 3; //was 5
//static int NUM_PRINTERS = 1;
static bool SAVE_STATE=false;

static int FILE_COUNTER=0;

static dGeomID floorplane;

static Network* cppn_network;
static bool first;

//static char startgenes_fn[100]="bipedstartgenes";
//static char substrate_fn[100] = "substrategenes";
static int novelty_function = 0;

//#define RIGIDITY 

dReal planes[]= // planes for a cube
  {
    1.0f ,0.0f ,0.0f ,0.05f,
    0.0f ,1.0f ,0.0f ,0.05f,
    0.0f ,0.0f ,1.0f ,0.05f,
    -1.0f,0.0f ,0.0f ,0.05f,
    0.0f ,-1.0f,0.0f ,0.05f,
    0.0f ,0.0f ,-1.0f,0.05f
    /*
    1.0f ,0.0f ,0.0f ,2.0f,
    0.0f ,1.0f ,0.0f ,1.0f,
    0.0f ,0.0f ,1.0f ,1.0f,
    0.0f ,0.0f ,-1.0f,1.0f,
    0.0f ,-1.0f,0.0f ,1.0f,
    -1.0f,0.0f ,0.0f ,0.0f
    */
  };
const unsigned int planecount=5; //was 6

dReal points1[]= // points for a cube
  {
    +0.25f,-0.25f,0.05f,  //  point 0
    -0.25f,+0.25f,0.05f, //  point 1
    +0.25f,+0.25f,0.05f, //  point 2

    0.25f,-0.25f,-0.05f,  //  point 3
    -0.25f,0.25f,-0.05f, //  point 4
    +0.25f,0.25f,-0.05f, //  point 5

  };

dReal points2[]= // points for a cube
  {
    -0.25f,-0.25f,0.05f,  //  point 0
    -0.25f,+0.25f,0.05f, //  point 1
    +0.25f,+0.25f,0.05f, //  point 2

    -0.25f,-0.25f,-0.05f,  //  point 3
    -0.25f,0.25f,-0.05f, //  point 4
    +0.25f,0.25f,-0.05f, //  point 5

  };

const unsigned int pointcount=6; //was 8
unsigned int polygons[] = //Polygons for a cube (6 squares)
  {
    3,2,1,0, // positive X Side 0
    3,3,4,5, // negative Y Side 4
    4,1,4,3,0, // positive Y Side 1    0, 3, 4, 1
    4,2,5,4,1, // positive Z Side 2   1, 4, 5, 2
    4,0,3,5,2, // negative X Side 3
  
   // 4,5,4,6,7, // negative Z Side 5
  };

#define NF_COG 14
#define NF_COGSAMP 15
#define NF_LEG 16
#define NF_LEGSAMP 17
#define NF_SS 18

#define NF_FITCU 0
#define NF_FITCUSAMP 1
#define NF_COGCU 2
#define NF_COGCUSAMP 3
#define NF_LEGCU 4
#define NF_LEGCUSAMP 5
#define NF_COGSQ 6
#define NF_COGSAMPSQ 7
#define NF_LEGSQ 8
#define NF_LEGSAMPSQ 9
#define NF_FITSQ 10
#define NF_FITSQSAMP 11
#define NF_FIT 12
#define NF_FITSAMP 13

//*********************
//#define ST_FS_LEO 0 
//#define ST_ES 1
//#define ST_ES_LEO 2

//#define DEBUG_OUTPUT

using namespace std;

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <ode/ode.h>

#ifdef _MSC_VER
#pragma warning(disable:4244 4305)  // for VC++, no precision loss complaints
#endif

#ifdef GRAPHICS
#include <drawstuff/drawstuff.h>
#include "texturepath.h"
#endif

class WalkerDomain;
static WalkerDomain *domain;
//static SGA::Population *pop;
//static SGA::GA *ga;
static dReal params[]={ 0.2, 3.0, 0.1,     0.15, 0.0, 0.1 ,   0.2,0.3,0.35};

static bool INItIAL_COLLESION=false;
static int COLLISION_COUNTER=0;

class SineController;
class DummyController;
class Controller;

static Controller* controller;

//NEAT + NS stuff
inline float dist(float x1, float y1, float x2, float y2)
{
	float xd = x1-x2;
	float yd = y1-y2;
	return xd*xd+yd*yd;
}

static	void calculate_delta(dVector3 v1, dVector3 v2, dVector3 o)
{
	for(int x=0;x<3;x++)
		o[x]=v2[x]-v1[x];
}

static	void calculate_power(dVector3 v,int pow)
{
	for(int x=0;x<3;x++)
	{
		float temp=v[x];
		bool sign=false;
		if(temp<0.0)
			sign=true;
		for(int k=1;k<pow;k++)
			v[x]*=temp;
		if(sign)
			v[x]=(-v[x]);
	}
}
float feet_distance2(vector<float>& t1, vector<float>& x1, vector<float>& y1, vector<float>& t2, vector<float> &x2, vector<float> &y2)
{
	float td=0.0;
	float dummy_distance=0.05;
	int c1=0;
	int c2=0;
	int size1 = x1.size();
	int size2 = x2.size();
	int loopsize=0;
	int bigsize=0;
	float  step_factor=1.0;
	float inc_factor=1.0;

	//cout << x1.size() << " " << y1.size() << " " << x2.size() << " " << y2.size() << endl;
	if(size1<size2)
	{
		loopsize=size1;
		bigsize=size2;	
	}
	else
	{
		loopsize=size2;
		bigsize=size1;
	}
	for(int x=0;x<loopsize;x++)
	{
		//	cout << x << " " <<  y1.size() << " " << y2.size() << endl;
		td+=step_factor*dist(x1[x],y1[x],x2[x],y2[x]);
		//step_factor*=inc_factor;
	}

	for(int x=loopsize;x<bigsize;x++)
	{
		//	cout << x << " " << loopsize << " " << bigsize << endl;
		if(size1<size2)
			td+=step_factor*dist(x2[x],y2[x],0.0,0.0);
		else
			td+=step_factor*dist(x1[x],y1[x],0.0,0.0);
		//step_factor*=inc_factor;
	}

	return td;
}

float feet_distance(vector<float>& t1, vector<float>& x1, vector<float>& y1, vector<float>& t2, vector<float> &x2, vector<float> &y2)
{
	float td=0.0;
	int c1=0;
	int c2=0;
	int size1 = x1.size();
	int size2 = x2.size();

	int time=0;

	float  step_factor=1.0;
	int step_cnt=0;
	float dummy_distance=0.05;
	float inc_factor=1.00;

	while(c1<size1 && c2<size2)
	{
		if(t1[c1+1]<t2[c2+1])
		{
			td+=step_factor*dist(x1[c1],y1[c1],x2[c2],y2[c2])*(t1[c1+1]-time);
			time=t1[c1+1];
			c1++;
		}
		else if (t2[c2+1]<t1[c1+1])
		{
			td+=step_factor*dist(x1[c1],y1[c1],x2[c2],y2[c2])*(t2[c2+1]-time);
			time=t2[c2+1];
			c2++;
		}
		else	
		{
			td+=step_factor*dist(x1[c1],y1[c1],x2[c2],y2[c2])*(t2[c2+1]-time);
			time=t2[c2+1];
			c1++;
			c2++;
		}
		if(c1 > step_cnt)
		{
			step_cnt=c1;
			step_factor*=inc_factor;
		}
		else if(c2>step_cnt)
		{
			step_cnt=c2;
			step_factor*=inc_factor;
		}
	}

	vector<float>* t;
	vector<float>* x;
	vector<float>* y;
	int c;

	if (t1[t1.size()-1] < t2[t2.size()-1])
	{
		t=&t2;
		x=&x2;
		y=&y2;
		c=c2;
	}
	else
	{
		t=&t1;
		x=&x1;
		y=&y1;
		c=c1;
	}

	step_cnt=c;
	step_factor=pow(inc_factor,c);

	for(int k=0;k<t->size()-1;k++)
	{
		td+=step_factor*dummy_distance*((*t)[k+1]-time);
		time=(*t)[k+1];
		step_factor*=inc_factor;
	}

	if(td<0.0)
	{
		cout << "Wtf" << endl;
		return 0.001;
	}
	return td;	
}

//novelty metric for maze simulation
float walker_novelty_metric(noveltyitem* x,noveltyitem* y)
{
	float dist=0.0;
	if(novelty_function!=NF_SS)
	{

		int size = x->data[0].size();

		for(int k=0;k<size;k++)
		{
			float delta = x->data[0][k]-y->data[0][k];
			dist+=delta*delta;
		}

		return dist;	

	}
	else
	{
		float left_feet;
		float right_feet;
		left_feet = feet_distance2(x->data[0],x->data[1],x->data[2],y->data[0],y->data[1],y->data[2]);
		right_feet = feet_distance2(x->data[3],x->data[4],x->data[5],y->data[3],y->data[4],y->data[5]);

		return left_feet+right_feet;
	}
}

void biped_epoch(NEAT::Population *pop,bool novelty=false);
void biped_realtime_loop(NEAT::Population *pop,bool novelty=false);

NEAT::Population *biped_realtime(bool novelty=false);
noveltyitem* biped_evaluate(NEAT::Organism* org,data_record* data=NULL);

//globals
static NEAT::Population *neatpop;
static noveltyarchive archive(1.0,*walker_novelty_metric, false); //! was 1
static data_rec Record; //stores run information
static int indiv_counter=0;
static bool Novelty=false;
static char outdir[50]="";

static Genome* substrate_genome;

static dReal DENSITY=0.1; //was 0.5  
//static dReal MAX_PARTS=250;

static dReal P_CONSTANT= 9.0; //was 9.0
static dReal D_CONSTANT= 0.0;
static dReal FOOTFACTOR= 5.0;

static bool increase_leg_length = false;
static bool decrease_leg_length = false;

static bool bDisplay = false;
static bool bDoEvolution = true;
static bool bSlowDown = false;
#ifdef GRAPHICS
static bool bNoVis=false; //true;
#else
static bool bNoVis=true;
#endif
static bool bMoviePlay=false;


dReal scaleval(dReal v, dReal min, dReal max)
{
	if(v<min)
		return 0.0;
	if(v>max)
		return 1.0;
	return (v-min)/(max-min);
}

class Biped: public Creature
{
public:
	dWorldID _worldID;
	dSpaceID _spaceID;
	int partNumber;
	int id1, id2;
	int timer;
	dVector3 orig_com;
	dVector3 orig_left;

	float tmpTimer;

	dBodyID bodID;
	int boxID;
	vector<dBodyID> tempbody;
	//dBodyID tempbody;
	int markerID;
	int nozzleID;
	dJointID tempjoint;
	bool first;
	vector< dVector3* > oldPosition;
	double energyConsumption;


	vector< vector<dGeomID> > machines;
	vector<int> bends;
	vector<double> latitudeAngles;
	vector<double> longitudeAngles;
	dQuaternion rightUp;
	dQuaternion rightDown;
	dQuaternion leftUp;
	dQuaternion leftDown;

	int anchorJointID;

	vector<dVector3*> jointPosition;
	dGeomID convexID;

	dGeomID closestGeom;

	vector<int> jointIDs;
	//, joint2;
	dJointID jointID;

	~Biped()
	{
	}

		
	
	static void initialNearCallback (void *data, dGeomID o1, dGeomID o2)
	{
		dBodyID b1,b2;
		//convexID=NULL;

		b1 = dGeomGetBody(o1);
		b2 = dGeomGetBody(o2);
		
		if(o1 == floorplane || o2 == floorplane)
			return;

		//TODO what we should really check is if o1 and o2 are next to each other
		//Different than GECCO 
		//if (b1 && b2 && dAreConnectedExcluding (b1,b2,dJointTypeContact)) return;
	//we don't know if they objects are connected yet
		//	if (b1 && b2 && dAreConnected(b1,b2)) return;
		//b1.

		const dReal* pos1 = dGeomGetOffsetPosition(o1);
		const dReal*  pos2 = dGeomGetOffsetPosition(o2);
		float dist = sqrt ( pow(pos1[0]-pos2[0],2) + pow(pos1[1]-pos2[1],2) +pow(pos1[2]-pos2[2],2));
		
		if (dist>=0.45) return;

	//	std::cout<<"COLLIDED WHILE BUILDING\t"<<pos1[0]<<"\t"<<pos2[0]<<"\t"<<dist<<"\n";
		INItIAL_COLLESION=true;
		return;

		const int N = 32;
		dContact contact[N];
		int n = dCollide (o1,o2,N,&(contact[0].geom),sizeof(dContact));
	//cout<<n<<"\n";
		if (n > 0) 
		{
			INItIAL_COLLESION=true;
		}

	}


	Biped(bool logging=false,bool movie=false):Creature(logging,movie) {
		partNumber=0;
		tmpTimer = 0.0f;
		first = true;
	//	convexID = new dGeomID(-1);
		
		energyConsumption = 0;

		//hack
		for (int i=0; i<30; i++)///NEAT::NumMachines; i++)
		{
			tempbody.push_back(0);
			machines.push_back(vector<dGeomID>());
			jointIDs.push_back(0);
		}
		//bends=new vector<int>();

		float degreeMod=1.0;

		//right up
			dQuaternion Q2;
			dQFromAxisAndAngle(Q2, 0.0, 0, 1,  0.78);//
			dQuaternion Q3;
			dQFromAxisAndAngle(Q3, 1.0, 0, 0,  0.78*degreeMod);// 
			dQuaternion Q4;
			dQuaternion Q5;
			dQMultiply0(Q5,Q3,Q2);
			dQFromAxisAndAngle(Q2, 0.0, 0, 1,  -0.78);//
			dQMultiply0(rightUp,Q2,Q5);
					
			//right down	
			dQFromAxisAndAngle(Q2, 0.0, 0, 1,  0.78);//
			dQFromAxisAndAngle(Q3, 1.0, 0, 0,  -0.78*degreeMod);//
			dQMultiply0(Q5,Q3,Q2);
			dQFromAxisAndAngle(Q2, 0.0, 0, 1,  -0.78);//
			dQMultiply0(rightDown,Q2,Q5);
					
			//left up 
			dQFromAxisAndAngle(Q2, 0.0, 0, 1,  -0.78);//
			dQFromAxisAndAngle(Q3, 1.0, 0, 0,  0.78*degreeMod);//
			dQMultiply0(Q5,Q3,Q2);
			dQFromAxisAndAngle(Q2, 0.0, 0, 1,  +0.78);//
			dQMultiply0(leftUp,Q2,Q5);
						

			//left down
			dQFromAxisAndAngle(Q2, 0.0, 0, 1,  -0.78);//
			dQFromAxisAndAngle(Q3, 1.0, 0, 0,  -0.78*degreeMod);//
			dQMultiply0(Q5,Q3,Q2);
			dQFromAxisAndAngle(Q2, 0.0, 0, 1,  +0.78);//
			dQMultiply0(leftDown,Q2,Q5);
					
	} 


	virtual dReal fitness()
	{
/*
#ifdef GRAPHICS
		return 0;
#endif
*/
		dVector3 centermass;

		CenterOfMass(centermass);
#ifdef GRAPHICS
		dMatrix3 RI;
		dRSetIdentity (RI);
	//	dVector3 v = {centermass[0], centermass[1], centermass[2]};
		const dReal ss2[3] = {5.02,5.02,5.02};

		//dsDrawBox (centermass, RI, ss2);
#endif

		
		//return sqrt( pow(dGeomGetPosition(closestGeom)[0]-startpos[0],2)+pow(dGeomGetPosition(closestGeom)[1]-startpos[1],2));
		if (WATER)
		{
			return sqrt( pow(centermass[0]-startpos[0],2)+pow(centermass[1]-startpos[1],2)+ pow(centermass[2]-startpos[2],2));
		}

//		if (energyConsumption==0)
//			return 0.0;

		//!double ff=sqrt( pow(centermass[0]-startpos[0],2)+pow(centermass[1]-startpos[1],2))/ (1.0+energyConsumption/5000.0);
			
		double ff = sqrt( pow(centermass[0]-startpos[0],2)+pow(centermass[1]-startpos[1],2)) - sqrt( pow(centermass[2]-startpos[2],2));
			
		//
		//if (INItIAL_COLLESION)
		//{
		//	ff/=10.0f;
		//	INItIAL_COLLESION=false;
		//}
		return ff;

		//Space filling fitness
		dVector3 center = {1.0, 0.0, 30.0};// {printerPos, 0.0, z_level};

		float fitness=0.0;
		float dx = 4.0;
		for (float xp=center[0]-dx; xp<center[0]+dx; xp+=0.5)
				for (float yp=center[1]-dx; yp<center[1]+dx; yp+=0.5)
						for (float zp=center[2]-dx; zp<center[2]+dx; zp+=0.5)
						{
				
								//bool close=false;
							for (int i=4; i<bodies.size(); i++)
							{
							//	if (i!=id1 && i!=boxID && i!=nozzleID)
							//	{
									const dReal *bp = dBodyGetPosition(bodies[i]);
									if (sqrt( pow(bp[0]-xp,2)+pow(bp[1]-yp,2)+pow(bp[2]-zp,2))<0.5)
									{
										fitness++;
										break;
									}
								//}
							}
							//cout<<"Test";

						}
	//	cout<<"\nFITNESS "<<fitness<<"\n";

		return fitness;
	}

	//~Biped()
	//{
	//	if(log)
	//		delete logfile;
	//}

	virtual bool abort()
	{
		//return false;
		//dVector3 pos;
		//this->CenterOfMass(pos);
		////cout<<startpos[2]<<"\n";
		//if (pos[2]>15 || startpos[2]>3.0)
		//{
		//	return true;
		//}

		//if (timer==100)
		//{
		//	dVector3 centermass;
		//	CenterOfMass(centermass);
		//   
		//	//cout<<sqrt( pow(centermass[0]-startpos[0],2)+pow(centermass[1]-startpos[1],2))<<"\n";
		//	if ( sqrt( pow(centermass[0]-startpos[0],2)+pow(centermass[1]-startpos[1],2))<0.5)
		//   {
		//	   return true;
		//   }
		//}

		if (INItIAL_COLLESION) 
			{
				INItIAL_COLLESION=false;
				return true;
		}

		return false;
	}


	//http://www.supermanoeuvre.com/blog/?p=671

	virtual void Create(dWorldID worldi, dSpaceID spacei, dVector3 posi,Controller* cont)
	{
		//SCALE_FACTOR = 1;
		_worldID = worldi;
		_spaceID = spacei;
		timer=0;
		dVector3 xAxis={1.0,0.0,0.0};
		dVector3 nxAxis={-1.0,0.0,0.0};
		dVector3 yAxis={0.0,1.0,0.0};
		dVector3 zAxis={0.0,0.0,1.0};

		Creature::Create(worldi,spacei,posi,cont);
		//dVector3 upper_pos = {1, 1, 1};
		//int upperleg = add_cylinder(1,DENSITY,0.5,0.5,upper_pos);

	}


	void resetWorld()
	{
		Destroy();

		if (CPPN_DET_NUM_PARTS)
			NEAT::NumParts=-1;

		bodies.clear();
		geoms.clear();
		joints.clear();

		//Destroy();
		Create(_worldID,_spaceID,pos,0);
		partNumber=0;
		for (int i=0; i<NEAT::NumMachines; i++)
			machines[i].clear();
	}

	int old_sec;


	void checkForNewGenome()
	{
		struct tm* clock;				// create a time structure

		struct stat attrib;			// create a file attribute structure

		char file_name[50];
		sprintf(file_name,"genome%d.txt",NEAT::REQUERY);

		//ostringstream file_name;
		//file_name.str("");
		//file_name<<"genome"<<NEAT::REQUERY<<".txt";


		if (stat(file_name, &attrib)!=0)
		{
			return;		// get the attributes of afile.txt
		}

		clock = gmtime(&(attrib.st_mtime));
		if (first) //first init
		{
			old_sec =clock->tm_sec;
			
		//	srand( (unsigned)time( NULL ) );

			//float xyz[3] = {-1.6,0.4,31.02};
			//float hpr[3] = {-25,-27, 0};
			//dsSetViewpoint(xyz, hpr);
			static float xyz[3] = {-35.6,9.9,38.3};
			static float hpr[3] = {-29,-23.5,0};
			#ifdef GRAPHICS
					dsSetViewpoint (xyz,hpr);
			#endif
		}
		first=false;
		//	cout<<clock->tm_sec<<"\n";
		if (old_sec!=clock->tm_sec) //TODO sec could be the same
		{
			int id=0;
			Genome* start_genome;

			//cout<<"File changed";
			ifstream iFile(file_name, ios::in);
			start_genome = new Genome(id, iFile);
			iFile.close();

			// Create CPPN from start genome
			//old_cppn_network = cppn_network;
			cppn_network = start_genome->genesis(id);
			dVector3 pos={0.0,0.0,0.0};

			resetWorld();
		}
		old_sec = clock->tm_sec;

		//WIN32_FILE_ATTRIBUTE_DATA fileAttrData = {0};
		//GetFileAttributesExW(L"C:\\abc.txt, GetFileExInfoStandard, &fileAttrData);

	}

	//http://en.wikipedia.org/wiki/STL_(file_format)
	//http://www.sv.vt.edu/classes/vrml/NCSA_VRML_Tutorial/examples/Walls.wrl.txt
	//http://www.ode.org/old_list_archives/2007-October/022922.html ONLY COLLIDE WITH THE GROUND
	//  dBodyID b1 = dGeomGetBody(o1);
  //dBodyID b2 = dGeomGetBody(o2);
 // if (b1 && b2 && dAreConnectedExcluding (b1,b2,dJointTypeContact)) return;
	
	//http://x3dgraphics.com/examples/X3dForWebAuthors/Chapter02-GeometryPrimitives/
	//http://www.adrian.zentner.name/content/projects/xml/x3d/rotation/index.html
	//TRY: Novelty assistant interactive evolution
	dBodyID createShape(int machineNumber, float printerPos, dBodyID tempbody, bool& finished, int maxParts=NEAT::NumParts, bool showPrinter=true)
	{
		double z_level = 30;
	
		int currentCount = machines[machineNumber].size();
		//cout<<currentCount<<"\n";

		while (machines[machineNumber].size()< (currentCount+maxParts) )
		{
		
			if (machines[machineNumber].size() < (currentCount+maxParts))//&& (timer % 5==0) )
			{
				//TODO don't creat them again every time
				dVector3 xAxis={1.0,0.0,0.0};
				dVector3 zAxis={0.0,0.0,1.0};

				double circ = 0.05; //was0.05
				double segment_size = 0.25; // 0.3 0.25
				if (machines[machineNumber].size()==0)
				{

					latitudeAngles.clear();
					longitudeAngles.clear();
					bends.clear();
					dVector3 posBox = {printerPos/2.5,-50.0, 15.0};//was -20.0

					//showPrinter=false;
					if (showPrinter)
					{
						boxID = add_box(100, 2.0, 15.0, 30.0, posBox); 

					//markerID = add_box(0.0001, 0.3,0.3,0.3,posBox);


						dVector3 orientation={0.0, 1.0,0.0};
						dVector3 pos2 = {printerPos/2.5, -20.5, z_level}; //was -10.5
						//was 0.3
						nozzleID = add_cylinder(orientation,DENSITY,40.0,1.0,pos2);  //was 20.0
						add_fixed(nozzleID, boxID);
					}

				}
				else
				{
					//	dJointDestroy(tempjoint);
				}
			//	dGeomID maingeom = dCreateGeomTransform(space);
			
				partNumber++;
				dVector3 newpos = {0.0,0.0,0.0}; //was -1.0
				dVector3 orientation={0.0,1.0,0.0}; //0.0, 1.0, 0.0


				dGeomID maingeom;


				double coin = (double)(rand()/(RAND_MAX+1.0));

				double zv=0.78; //0.78
				double yv = 0.78;

				double inputs[CPPN_INPUTS];

				//GECCO SETTINGS
				inputs[0]=  (machines[machineNumber].size() - (NEAT::NumParts-1) / 2.0) / ( (NEAT::NumParts-1) / 2.0) ; 
				inputs[1]= abs(inputs[0]); //BIAS // (machines[machineNumber].size() - (NEAT::NumParts-1) / 2.0) / ( (NEAT::NumParts-1) / 2.0) ;
				
				//Setting for bipedal video
				//inputs[0]=  (machines[machineNumber].size() - (NEAT::NumParts-1) / 2.0) / ( (NEAT::NumParts-1) / 2.0) ;
				//inputs[0]= abs(inputs[0]);// (machines[machineNumber].size() - (NEAT::NumParts-1) / 2.0) / ( (NEAT::NumParts-1) / 2.0) ;
				//inputs[1]=0;
				//inputs[2] = 1.0;//bias

				//SYMMETRY
			//	inputs[1] = abs(inputs[0]);
			//	inputs[0] = 0;

				//if (inputs[0]<0.) inputs[0]*=-1.0;

				
			//	if (CPPN_INPUTS>3)
			//		inputs[3] = printerPos;

				//cout<<inputs[0]<<"\t";

			
				cppn_network->flush();
				cppn_network->load_sensors(inputs);
				//cout<<cppn_network->nodecount()<<"\n";
				//CPPN_QUERIES
				//!cppn_network->nodecount();  //was 6 (1.22.13)
				for(int i = 0; i < CPPN_QUERIES; i++) //cppn_network->nodecount(); i++) //was 15fUpda
					cppn_network->activate();

				int index=-1;
				double max = -10000;

				zv = 0.0;
				yv = 0.0;
				float dif=0.0f;
				float dif2=0.0f;
				float dif3=0.0;
				float shouldBend = 0.0f;

				if ( cppn_network->outputs.size() >0)
				{
					dif=cppn_network->outputs[0]->activation*2.0-1.0;//-cppn_network->outputs[1]->activation;
					dif2=cppn_network->outputs[1]->activation*2.0-1.0;//-cppn_network->outputs[3]->activation;
				//	dif2=cppn_network->outputs[1]->activation;
					
					dif3=cppn_network->outputs[1]->activation*2.0-1.0;
				
				//	cout<<dif<<"\n";
					//shouldBend = cppn_network->outputs[2]->activation*2.0-1.0;
				}

				
				int bend=4; //No bend

				//if (shouldBend > 0.4)
				if (abs(dif)>abs(dif2))
				{
					if (abs(dif)>0.5)
					{
						if (dif<0)
						{
							zv=0.78;
							bend = 0;

						}
						else
						{
							zv=-0.78;
							bend = 2;
						}
					}								
				}
				else if (abs(dif2)>abs(dif3))
				{

					if (abs(dif2)>0.5)
					{
						if (dif2<0)
						{
							yv=0.78;
							bend = 1; //up
						}
						else
						{
							yv=-0.78;
							bend = 3; //down
						}
					}
				}
				else
				{
					if (abs(dif2)>0.5)
					{
						if (dif2<0)
						{
							yv=0.78;
							bend = 5; //up
						}
						else
						{
							yv=-0.78;
							bend = 6; //down
						}
					}
				}

				//if (dif>0)
				//{
				//	if (dif<0.5)
				//		bend = 0;
				//	else
				//		bend = 1;
				//}
				//else
				//{
				//	if (dif>-0.5)
				//		bend = 2;
				//	else
				//		bend = 3;
				//}

				//if (dif>0.5)
				//{
				//	bend = 3;
				//}
				//else bend = 4; //no bend

				dQuaternion* Q;

				//bend = manualBend[r];
				//bend=4;

			//	bend = (int)(machines[machineNumber].size()/50);
			//	bend=36;
				//machines[machineNumber].size()==1 || 

				//if (machines[machineNumber].size()>80 && machines[machineNumber].size()<90)
				//	bend = 4;

				//bend=1;

				//Part before and after motor should be straight
				if (machines[machineNumber].size()>currentCount+maxParts-5 || (machines[machineNumber].size()<=currentCount+2) )
				{
					//bend = 4;
				}
				//	bend=(int)(rand()%5);

				dGeomID convex0=NULL;

				//if (machines[machineNumber].size()%6==0)
				//{
				//	bend = 1;
				//}
				//else
				//{
				//	bend = 3;
				//}

				//bend = 3;

				if (bend==3||bend==2||bend==4||bend==5||bend==6)
				{
					convex0= dCreateConvex (space,  planes,  planecount,  points2,  pointcount, polygons);
					dGeomSetData(convex0, &points2);
				}
				else
				{
					convex0= dCreateConvex (space,  planes,  planecount,  points1,  pointcount, polygons);
					dGeomSetData(convex0, &points1);
				}

				geoms.push_back(convex0);
				dGeomID convex1;
				if (bend==3||bend==2||bend==4||bend==5||bend==6)
				{
					convex1= dCreateConvex (space,  planes,  planecount,  points2,  pointcount, polygons);
					dGeomSetData(convex1, &points2);
				}
				else
				{
					convex1= dCreateConvex (space,  planes,  planecount,  points1,  pointcount, polygons);
					dGeomSetData(convex1, &points1);
				}

				//dCreateConvex (space,  planes,  planecount,  points,  pointcount, polygons);
				geoms.push_back(convex1);

				dGeomSetBody(convex0,tempbody);
				dGeomSetBody(convex1,tempbody);
				
				dGeomSetOffsetPosition(convex0, newpos[0], newpos[1], newpos[2]);
				dGeomSetOffsetPosition(convex1, newpos[0], newpos[1], newpos[2]);

				dQuaternion Q2;
				dQFromAxisAndAngle(Q2, 0.0, 0, 1,  3.1415);//
				dGeomSetOffsetQuaternion(convex1, Q2);
			

				if (bend==0) //right up
				{
					Q = &rightUp;
					dGeomSetOffsetQuaternion(convex0, rightUp);
				}
				else if (bend==1) //right down
				{
					Q=&rightDown;
					dGeomSetOffsetQuaternion(convex0, rightDown);
				}
				else if (bend==2) //left up 
				{
					Q=&leftUp;
					dGeomSetOffsetQuaternion(convex0, leftUp);
				}

				else if (bend==3) //left down
				{
					Q=&leftDown;
					dGeomSetOffsetQuaternion(convex0, leftDown);
				}
				else if (bend==5)
				{
					dQuaternion m;
					Q=&m;
					dQFromAxisAndAngle( m, 1, 0, 0,  0.78);
					dQuaternion qold;
					dQuaternion Q2;
					dGeomGetOffsetQuaternion(convex0, qold);
					dQMultiply0(Q2,(*Q),qold);
					dGeomSetOffsetQuaternion(convex0, Q2);

					dGeomGetOffsetQuaternion(convex1, qold);					
					dQMultiply0(Q2,m,qold);
					dGeomSetOffsetQuaternion(convex1, Q2);

					
				}
				else if (bend==6)
				{
					dQuaternion m;
					Q=&m;
					dQFromAxisAndAngle( m, 1, 0, 0,  -0.78);
					dQuaternion qold;
					dQuaternion Q2;
					dGeomGetOffsetQuaternion(convex0, qold);
					dQMultiply0(Q2,(*Q),qold);
					dGeomSetOffsetQuaternion(convex0, Q2);

					dGeomGetOffsetQuaternion(convex1, qold);					
					dQMultiply0(Q2,m,qold);
					dGeomSetOffsetQuaternion(convex1, Q2);

					
				}

			
				machines[machineNumber].push_back(convex0);
				machines[machineNumber].push_back(convex1);

				bends.push_back(bend);


		
				//MOVE OLD PIECES TO THEIR NEW POSITIONS
				for (int i=machines[machineNumber].size()-3; i>=0; i--) //was -2
				{

					dQuaternion rotationQ;
					dQuaternion newQ2;
					dQuaternion newQ3;
					
					//dQuaternion Q;
					dQuaternion newQ;

					dGeomID current = machines[machineNumber][i];

					dGeomID previous= machines[machineNumber][i+1];

					dQuaternion oldQ;
					dGeomGetOffsetQuaternion(current, oldQ);


					if (bend!=4)
					{

						dQMultiply0(newQ2,  (*Q), oldQ);//R3);
						dNormalize4(newQ2);
						dGeomSetOffsetQuaternion(current,newQ2);
					
					}

					dVector3 *oldpos2 = new dVector3[3];
					(*oldpos2)[0] = dGeomGetOffsetPosition(previous)[0];
					(*oldpos2)[1] = dGeomGetOffsetPosition(previous)[1];
					(*oldpos2)[2] = dGeomGetOffsetPosition(previous)[2];


					dVector3 newpos;
					
					dVector3 beginning;

					if ( (i+1)%2==0)
					{ //{0.0,0.25, 0.00};
						//i=1

							dVector3 lengthVector = {0.00,0.25, 0.00};//{0.0,0.0,segment_size/2.0}; //{0.5,0.5,0.5};
							dMULTIPLY0_331(newpos, dGeomGetOffsetRotation(previous), lengthVector); //To world coordinates

					}
					else
					{

							dVector3 lengthVector = {0.00,0.25, 0.0};//{0.0,0.0,segment_size/2.0}; //{0.5,0.5,0.5};
							dMULTIPLY0_331(newpos, dGeomGetOffsetRotation(previous), lengthVector); //To world coordinates


					}

					if ( (i)%2==0)
					{

							dVector3 lengthVector = {0.00,0.0, 0.0};//{0.0,0.0,segment_size/2.0}; //{0.5,0.5,0.5};
							dMULTIPLY0_331(beginning, dGeomGetOffsetRotation(current),lengthVector); //Translate to world coordinates

					}
					else
					{ //0.0,-0.25, 0.0} {0.00,-0.25, 0.0};

							dVector3 lengthVector = {0.00,-0.25, 0.0};//{0.0,0.0,segment_size/2.0}; //{0.5,0.5,0.5};
							dMULTIPLY0_331(beginning, dGeomGetOffsetRotation(current),lengthVector); //Translate to world coordinates				
						
					}

					if (i%2==0)
					{
						dGeomSetOffsetPosition(current,dGeomGetOffsetPosition(previous)[0]+ beginning[0],dGeomGetOffsetPosition(previous)[1]+ beginning[1],dGeomGetOffsetPosition(previous)[2]+ beginning[2]);
					}
					else
					{
						dGeomSetOffsetPosition(current, (*oldpos2)[0]+newpos[0]+beginning[0], (*oldpos2)[1]+newpos[1]+beginning[1],(*oldpos2)[2]+newpos[2]+beginning[2]);
					}

					delete[] oldpos2;
				}
				

				
			//	dGeomID sphereID=dCreateSphere(space, 0.1);
			//	geoms.push_back(machines[machines.size()-1]);

				//cout<<"-------\n\n";
					
				//Check if there was a crash with the other machines

				//INItIAL_COLLESION=false;
				//dSpaceCollide (space,0,&initialNearCallback); 
				//if (INItIAL_COLLESION)
				//{
				//	cout<<"COLLISION\n";
				//	COLLISION_COUNTER++;
				//}
				//INItIAL_COLLESION=false;
		if (false)
		{
				INItIAL_COLLESION=false;
				dSpaceCollide (space,0,&initialNearCallback); //dWorldStep
			//INItIAL_COLLESION=true;

				if (INItIAL_COLLESION && bend!=4) //undue bend
				{
					dQuaternion qold;
					dGeomSetOffsetQuaternion(maingeom, qold);
					int counter=0;
					for (int i=machines[machineNumber].size()-2; i>=0; i--)
					{
			//!			dGeomSetOffsetQuaternion(machines[machineNumber][i],oldQuat[counter]);
			//!			dGeomSetOffsetPosition(machines[machineNumber][i],(*oldPos[counter])[0],(*oldPos[counter])[1] ,(*oldPos[counter])[2]);
						//cout<<(*oldPos[counter])[0]<<"\t";
						//cout<<(oldQuat[counter])[0]<<"\t";
						counter++;
					}
					//TODO delete elements in lists
					//cout<<"\n";
				}
				//for (int i=0; i<oldPos.size(); i++)
				//{
				//	delete[] oldPos[i];
				//	delete[] oldQuat[i];
				//}
		}
						//if (false && partNumber==NEAT::NumParts)
						//{
						//	int idx = geoms.size()-1;
						//	dVector3 lengthVector = {0.0,0.0,segment_size/2.0}; //{0.5,0.5,0.5};

						//	dVector3 newpos2;
						//	dMULTIPLY0_331(newpos2, dGeomGetOffsetRotation(geoms[idx]), lengthVector); //FOR VISUAL
						//	dVector3 oldpos3 = {dGeomGetOffsetPosition(geoms[idx])[0]+dBodyGetPosition(tempbody)[0], 
						//					dGeomGetOffsetPosition(geoms[idx])[1]+dBodyGetPosition(tempbody)[1],
						//					dGeomGetOffsetPosition(geoms[idx])[2]+dBodyGetPosition(tempbody)[2]}; //{2.0, 2.0, 2.0};
					
						//

						//						//FOR VISUAL
						//	dBodySetPosition(bodies[markerID], oldpos3[0]+newpos2[0],oldpos3[1]+newpos2[1],oldpos3[2]+newpos2[2]);
						//}


			}
			if (STEP_MODE) break;
		}
		//	INItIAL_COLLESION=false;
		//dSpaceCollide (space,0,&initialNearCallback);
		//if (INItIAL_COLLESION)
		//{
		//	cout<<"STILL COLLIDING";
		//}

		//	for (int i=0; i<bodies.size(); i++)
		//		dBodySetPosition (bodies[i], dBodyGetPosition(bodies[i])[0], dBodyGetPosition(bodies[i])[1]-0.01f, dBodyGetPosition(bodies[i])[2]);

		//if (machines[machineNumber].size()==NEAT::NumParts && !finished)
		//{
			finished=true;
		//adjustMass(machineNumber, tempbody, 0, maxParts-1, z_level, printerPos);
		//adjustMass(0, tempbody[m], start,start+ NEAT::NumParts/NEAT::NumMachines-1, ZLEVEL, 0.0);
		//}

		return tempbody ;

	}

	void adjustMass(int machineNumber, dBodyID tempbody, int start, int end, int z_level, float printerPos)
	{
		//return; //!
		//if (machines[machineNumber].size()==NEAT::NumParts && !finished)
		//{
			//finished = true;
			dMass m2;
			dMass m;
			dMassSetZero (&m);

			
			//int mp = machines[machineNumber].size()-maxParts;
			//machines[machineNumber].size()-1, mp
			for (int i=end; i>=start; i--)
			{
				dGeomID current = machines[machineNumber][i];

				//			dMassSetBox (&m2, /* Create normal mass for a box */
				//	0.1, 
				//	0.4,0.1,0.5); //! replace last entry with segment size 


				dMassSetBox (&m2, /* Create normal mass for a box */
					0.25,
					0.25,0.1,0.5); //! replace last entry with segment size 

				//dMassSetBox (&m2, /* Create normal mass for a box */
				//	0.3,
				//	0.4,0.1,0.5); //! replace last entry with segment size 

				dMassRotate (&m2, dGeomGetOffsetRotation(current));
				dMassTranslate (&m2,dGeomGetOffsetPosition(current)[0], 
													dGeomGetOffsetPosition(current)[1],
													dGeomGetOffsetPosition(current)[2]);

			//	dMassTranslate (&m,-m.c[0],-m.c[1],-m.c[2]);
				// add to the total mass
				dMassAdd (&m,&m2);

				//dGeomSetBody(geoms[i],tempbody);
			}

			//			dMatrix3 RI;
			//dRSetIdentity (RI);
			//const dReal ss[3] = {1.02,1.02,1.02};
			//dsDrawBox (m.c, RI, ss);

			for (int i=end; i>=start; i--)
			{
				dGeomID current = machines[machineNumber][i];
				dGeomSetOffsetPosition (current,dGeomGetOffsetPosition(current)[0]-m.c[0],dGeomGetOffsetPosition(current)[1]-m.c[1],	dGeomGetOffsetPosition(current)[2]-m.c[2]);
			}
						
		//	cout<<m.c[0]<<"\t"<<m.c[1]<<"\t"<<m.c[2]<<"\n";

			//dGeomID current = machines[machineNumber][machines[machineNumber].size()-1];
		
			//cout<<dGeomGetOffsetPosition(current)[2]<<"\n";
			//dBodySetPosition(tempbody , printerPos/2.5+dGeomGetOffsetPosition(current)[0]+m.c[0], 
				//									dGeomGetOffsetPosition(current)[1]+m.c[1],
					//								z_level+dGeomGetOffsetPosition(current)[2]+m.c[2]);



			dBodySetPosition(tempbody , printerPos/2.5+m.c[0],0.0+m.c[1],z_level+m.c[2]);


			dMassTranslate (&m,-m.c[0],-m.c[1],-m.c[2]);
			
		//	dMassTranslate (&m,-0.0,-0.0,-0.0);
			//dBodySetPosition(tempbody,printerPos/2.5,0.0,z_level); 

			dBodySetMass (tempbody ,&m);
			partNumber++;			
			//dVector3 po = {0.0,0.0,30.0}; //{0.5,0.5,0.5};

			//markerID = add_box(0.1, 1.3,1.3,1.3,po);
			//cout<<"\nFITNESS "<<fitness() <<"\n";
		//}
	}

	//COMPOSITE OBJECTS
	//http://www.panda3d.org/forums/viewtopic.php?t=9782

	string NumberToString ( float Number )
	{
		stringstream ss;
		ss << Number;
		return ss.str();
	}

	//http://support.ponoko.com/entries/21531613-how-to-export-a-stl-design-file-using-blender-2-6

	void actuateJoint(int jointID, double desired_angle)
	{
		float current_angles=dJointGetHingeAngle(joints[jointID]); 
		dReal delta=desired_angle-current_angles;
		double p_term = P_CONSTANT* delta;
		double d_term = (-D_CONSTANT*delta);
		desired_angle=p_term+d_term;
		dJointSetHingeParam(joints[jointID],dParamVel,desired_angle); 
	}

	void saveAsSCAD()
	{
		int m=0;
		if (SHAPE_MODE == 2)
			m = 1;
		else
			m = NEAT::NumMachines;

		for (int i=0; i<m; i++)
		{
			const dReal *pos;
			const dReal *R;
			string modelInfo="union() {\n";
			int counter=0;
			for(int x=0;x<NEAT::NumParts;x++)
			{
				if (counter==0)
				{
					modelInfo+="union() {\n";
				}
				 pos =dGeomGetPosition (machines[i][x]);// dGeomGetOffsetPosition (machines[i][x]);
				 dQuaternion q;
				// dQuaternion q2;

				 //const dReal *R;
				 dMatrix3 R2;
				  dMatrix3 R_final;


				 R=  dGeomGetRotation(machines[i][x]);//dGeomGetOffsetRotation(machines[i][x]);
			
				//dMULTIPLY0_333(R_final, R,R2);
			
				dRtoQ(R, q);
				// const dReal *q = dBodyGetQuaternion (bodies[x]);
						
			//	cout<<q[0]<<"\t"<<q[1]<<"\t"<<q[2]<<"\t"<<q[3]<<"\n";

				 double angel = 2.0 * acos(q[0]);
				//q[1] = 1;

				 double ax = angel!=0.0 ? q[1]/sin(angel/2.0) : 0.0;
				double ay = angel!=0.0 ? q[2]/sin(angel/2.0) : 0.0;
				double az = angel!=0.0 ? q[3]/sin(angel/2.0) : 0.0;
		
				angel = angel * (180.0/3.1415926535);
			 
					  //) { ")+" "+NumberToString(pos[1])+" "+NumberToString(pos[2]-28.0)+"'>\n";
			    		//	modelInfo+="<Transform  rotation ='1 0 0 1.57'>\n";

				//139
				  modelInfo+="translate(v = ["+NumberToString(pos[0])+", "+NumberToString(pos[1])+", "+NumberToString(pos[2])+"])"+"{\n";
				   modelInfo+="rotate(a = "+NumberToString(angel)+", v = ["+NumberToString(ax)+", "+NumberToString(ay)+", "+NumberToString(az)+"]){\n";
			
			 modelInfo+="   cube(size=[0.5,0.2,0.6],center=true);\n";
				 // modelInfo+=" <Box size='0.4 0.1 0.5'/>";
				 modelInfo+="}}\n";
				// modelInfo+="</Transform>\n";
				 counter++;
				 if (counter==25)
				 {
					  modelInfo+="}\n\n"; //end union
					  counter=0;
				 }

			}
			modelInfo+="}\n";; 
			std::ofstream myFile;
			//string filename = "model";
			//filename+=i;filename+=".scad";
			char file_name[50];

			sprintf(file_name,"model%ld.scad",i);//NEAT::REQUERY);
		

			myFile.open(file_name);//"model.scad");
			myFile << modelInfo.c_str();
			myFile.close();	

		}
		cout<<"Saving complete\n";
	}

	/**
	void saveAsX3D()
	{
		const dReal *pos;
		const dReal *R;
		string modelInfo="\n";
		
		 modelInfo+="<!DOCTYPE X3D PUBLIC \"ISO//Web3D//DTD X3D 3.1//EN\" \"http://www.web3d.org/specifications/x3d-3.1.dtd\">\n";
		 modelInfo+="<X3D profile='Interchange' version='3.1'  xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' xsd:noNamespaceSchemaLocation =' http://www.web3d.org/specifications/x3d-3.1.xsd '>\n";
		 modelInfo+="<head>\n";
		 modelInfo+="</head>\n"; 
		 modelInfo+="<Scene>\n";
		 modelInfo+="<Background skyColor='1 1 1'/> \n";
		 modelInfo+="<Viewpoint description='Book View' orientation='-0.747 -0.624 -0.231 1.05' position='-1.81 3.12 2.59'/>\n";


		for(int x=4;x<bodies.size();x++)
		{
			 pos = dGeomGetPosition (geoms[x]);
			 dQuaternion q;
			// dQuaternion q2;

			 //const dReal *R;
			 dMatrix3 R2;
			  dMatrix3 R_final;

			 R=  dGeomGetRotation(geoms[x]);
			
			//dMULTIPLY0_333(R_final, R,R2);
			
			dRtoQ(R, q);

			
			 float angel = 2.0 * acos(q[0]);
			//q[1] = 1;

			 float ax = angel!=0.0 ? q[1]/sin(angel/2.0) : 0.0;
			float ay = angel!=0.0 ? q[2]/sin(angel/2.0) : 0.0;
			float az = angel!=0.0 ? q[3]/sin(angel/2.0) : 0.0;

			// modelInfo+="<Transform translation='"+NumberToString(pos[1])+" "+NumberToString(pos[2]-28.0)+" "+NumberToString(pos[0])+"'>\n";
		
			 
			  modelInfo+="<Transform translation='"+NumberToString(pos[0])+" "+NumberToString(pos[1])+" "+NumberToString(pos[2]-28.0)+"'>\n";
			    	//	modelInfo+="<Transform  rotation ='1 0 0 1.57'>\n";

			   modelInfo+="<Transform rotation='"+NumberToString(ax)+" "+NumberToString(ay)+" "+NumberToString(az)+" "+NumberToString(angel)+"'\>\n";
			
			//  <Transform  rotation ='1 0 0 1.57'>
			// modelInfo+="<Transform rotation='"+NumberToString(ax)+" "+NumberToString(ay)+" "+NumberToString(az)+" "+NumberToString(angel)+"' translation='"+NumberToString(pos[0])+" "+NumberToString(pos[1])+" "+NumberToString(pos[2]-28.0)+"'>\n";
			 modelInfo+="<Shape>\n"; //0.4,0.1
			 modelInfo+=" <Box size='0.5 0.2 0.6'/>";
			 // modelInfo+=" <Box size='0.4 0.1 0.5'/>";
			 modelInfo+="<Appearance>\n";
			 modelInfo+="<Material/>\n";
			 modelInfo+="</Appearance>\n";
			 modelInfo+="</Shape>\n";
			 modelInfo+="</Transform>\n";
			// modelInfo+="</Transform>\n";
			 modelInfo+="</Transform>\n";
		}
		modelInfo+="</Scene>\n";
		modelInfo+="</X3D>\n"; 
		std::ofstream myFile;
		myFile.open("model.x3d");
		myFile << modelInfo.c_str();
		myFile.close();	
		cout<<"Saving complete\n";
	}*/

	// 2 input
	// ode/test/test_boxstack.cpp
	virtual void Update(double timestep)
	{
		

		if (RESET)
		{
			RESET = false;
			resetWorld();
		}


		if (SAVE_STATE)
		{
			saveAsSCAD();
			//saveAsX3D();
			//saveAsSTL();
			SAVE_STATE= false;

		}

		#ifdef GRAPHICS
				if (timer%2==0)
				{
					checkForNewGenome();
				}
		#else
					//resetWorld();
				//	partNumber=0; //always rebuild if we have no graphics enabled
		#endif


		if (NEAT::NumParts==-1) //means that CPPN_DET_NUM_MACHINES ==true
		{
			//CPPN_OUTPUTS = 3;
			double inputs[CPPN_INPUTS];
			inputs[0]=0.0; //pos
			inputs[1]=0; //distance from center
			//inputs[2]=1;//bias

			cppn_network->flush();
			cppn_network->load_sensors(inputs);
			//cout<<cppn_network->nodecount()<<"\n";
			for(int i = 0; i < CPPN_QUERIES; i++) //cppn_network->nodecount(); i++) //was 15
				cppn_network->activate();

			//also change number of motors
			NEAT::NumParts = (int)( cppn_network->outputs[2]->activation*240.0+60.0);//-cppn_network->outputs[3]->activation;

			//If the number of parts is determined by the CPPN the give number of machines give the number for a 80 piece robot. Therefore if the robot is longer we have to scale it
			NEAT::NumMachines = (int)(SEGMENTS_BETWEEN_MOTORS *(NEAT::NumParts/60.0));

			//make sure numparts is divisible by nummachines
			NEAT::NumParts= ((int)(NEAT::NumParts/NEAT::NumMachines))*NEAT::NumMachines;
			
		//	cout<<cppn_network->outputs[2]->activation<<"\t"<<NEAT::NumParts<<"\n";
			//cout<<NEAT::NumParts<<"\t"<<NEAT::NumMachines<<"\n"; 
			//if (NEAT::NumParts<80)
			//{
			//	cout<<"Shouldn't happen";
			//}
		}

		if (machines[0].size()<NEAT::NumParts)// partNumber==0)
		{
			bool finished=false;

			if (NEAT::NumMachines==1)
			{
				///dBodyID tempbody;
				if (machines[0].size()==0)
				{
					tempbody[0]= dBodyCreate (world);
					dBodySetPosition(tempbody[0] ,0,0,30.0); //10 was zlevel
					dBodySetData (tempbody[0] ,(void*) 1);
				}
				createShape(0, 0.0,tempbody[0], finished);
			}
			else
			{
				if (SHAPE_MODE==1)
				{
					if (machines[0].size()==0)
					{

						for (int i=0; i<NEAT::NumMachines; i++)
						{
							tempbody[i]= dBodyCreate (world);
							dBodySetPosition(tempbody[i] ,0,0,30.0); //10 was zlevel
							dBodySetData (tempbody[i] ,(void*) 1);
						}
					}
					dBodyID id1 = createShape(0, -1.0, tempbody[0], finished); finished=false;
					dBodyID id2 = createShape(1, 0.0, tempbody[1], finished); finished=false;
					dBodyID id3 = createShape(2, 1.0, tempbody[2], finished);
					if (finished) //do something so we don't repeat this over and over again
					{
						//add_fixed(id1, id2);
						//add_fixed(id2, id3);
					}
				}
				else if (SHAPE_MODE==2)
				{
					dRandSetSeed(25);
					tmpTimer = 0;
					float ZLEVEL =30;


					if (machines[0].size()==0)
					{

						for (int i=0; i<NEAT::NumMachines; i++)
						{
							tempbody[i]= dBodyCreate (world);
							dBodySetPosition(tempbody[i] ,0,0, ZLEVEL); //10 was zlevel
							dBodySetData (tempbody[i] ,(void*) 1);

							bodies.push_back(tempbody[i]);
						}
					}
					if (machines[0].size()<NEAT::NumParts-1)
					{

						 finished=false;

						 for (int m=0; m<NEAT::NumMachines; m++)
						 {
							//dBodyID id1 = 
							 createShape(0, 0.0, tempbody[m], finished, NEAT::NumParts/NEAT::NumMachines,false); finished=false;
							//dBodyID id2 = 
							//createShape(0, 0.0, tempbody[1], finished, NEAT::NumParts/2.0); finished=false;
						 }


						
						//ADJUST MASS
						int start=0;
						for (int m=0; m<NEAT::NumMachines; m++)
						{

							adjustMass(0, tempbody[m], start,start+ NEAT::NumParts/NEAT::NumMachines-1, ZLEVEL, 0.0);
							start+= NEAT::NumParts/NEAT::NumMachines;
						}
						//adjustMass(0, tempbody[m], 0, NEAT::NumParts/2.0-1, 30.0, 0.0);
						//adjustMass(0, tempbody[1], NEAT::NumParts/2.0,  NEAT::NumParts-1, 30.0, 0.0);

						CenterOfMass(orig_com);
					//	CenterOfMass(curr_com);

						//Add joints
						start = NEAT::NumParts/NEAT::NumMachines;
						for (int m=0; m<NEAT::NumMachines-1; m++)
						{
							const dReal *pos = dGeomGetPosition(machines[0][start]); //NEAT::NumParts/3.0
							dVector3 lengthVector = {0.0,0.25,0.0}; //{0.5,0.5,0.5};
							dVector3 newpos;
							dMULTIPLY0_331(newpos, dGeomGetOffsetRotation(machines[0][start]), lengthVector); //To world coordinates
//							dVector3 pos2 = {newpos[0]+pos[0],newpos[1]+pos[1],newpos[2]+pos[2]};

							dVector3 anchor1 = {newpos[0]+pos[0],newpos[1]+pos[1],newpos[2]+pos[2]};//{pos[0],pos[1],pos[2]};//{,0,30}; //{0,0,30}
							

							dQuaternion q;
							dGeomGetOffsetQuaternion(machines[0][start], q);

							dVector3 yAxis;
							dVector3 axis = {0.0,1.0,0.0}; //was 0,0,1
					
							dMatrix3 R;
							//dRFromAxisAndAngle(R, 1, 0, 0, 1.52);
							dQtoR(q, R);
							dMULTIPLY0_331(yAxis, R,axis); //Translate to world coordinates
					
							//jointIDs[m] = add_hinge(tempbody[m],tempbody[m+1],anchor1,yAxis,-0.8,+0.8, 5.0);  //was 7.0 //was -1.4, 0.8
							jointIDs[m] = add_hinge(tempbody[m],tempbody[m+1],anchor1,yAxis,-1.57,+1.57, 5.0);  //was 7.0 //was -1.4, 0.8

							//1.57
							start+=NEAT::NumParts/NEAT::NumMachines;
						}

						oldPosition.clear();

						INItIAL_COLLESION=false;
					//	dSpaceCollide (space,0,&initialNearCallback); //dWorldStep
					//	if (INItIAL_COLLESION)
					//	{
					//		cout<<"collision\n";
					//	}

						for (int i=0; i<machines[0].size(); i++) //was -2
						{
							for (int e=0; e<machines[0].size(); e++) //was -2
							{
								if (i==e) continue;

								const dReal *pos1 = dGeomGetPosition(machines[0][i]);
								const dReal* pos2 = dGeomGetPosition(machines[0][e]);

								float dist = sqrt ( pow(pos1[0]-pos2[0],2) + pow(pos1[1]-pos2[1],2) +pow(pos1[2]-pos2[2],2));
								
								if (dist<0.45 && (dist!=0.00000)) //if it's exactly zero it's two triangles next to each other
								{
									//collision
									INItIAL_COLLESION=true;
									//cout<<"collision\t"<<dist<<"\n";
									break;
								}
							}
						}

						energyConsumption=0;
						for (int m=0; m<NEAT::NumMachines; m++)
						{
							for (int x=0; x<machines[m].size(); x++)
							{
								//dVector3 pos
								const dReal *pos = dGeomGetPosition(machines[m][x]);
								dVector3 *oldpos2 = new dVector3[3];
								(*oldpos2)[0]= pos[0];
								(*oldpos2)[1]= pos[1];
								(*oldpos2)[2]= pos[2];
								oldPosition.push_back(oldpos2);
							}
						}
						//CHECK IF BODIES INITIALLY COLLIDE WITH EACH OTHER - IF SO ABORT
	

						//initialNearCallback
			//Determine lowest z value
						int lowz=10000;
						//dVector3 startpos;
					
						for (int m=0; m<NEAT::NumMachines; m++)
						{
							for (int x=0; x<machines[m].size(); x++)
							{
								const dReal *pos = dGeomGetPosition(machines[m][x]); //NEAT::NumParts/3.0
								if (pos[2]<lowz)
								{
									lowz = pos[2];
									closestGeom = machines[m][x];

								}
							}

						}

						if (!WATER)
						{
							for (int m=0; m<bodies.size(); m++)
							{
								const dReal *pos = dBodyGetPosition(bodies[m]);
								dBodySetPosition(bodies[m],pos[0],pos[1],pos[2]-lowz);
								//const dReal *pos = dGeomGetPosition(machines[m][x]); //NEAT::NumParts/3.0
									//LMCDMA2010
					
							}
						}

						const dReal *pos2 = dGeomGetPosition(closestGeom); 

						startpos[0] = pos2[0];
						startpos[1] = pos2[1];
						startpos[2] = pos2[2];

						//cout<<"LOW: "<<lowz<<"\n";
						if (WATER)
						{
							dWorldSetGravity (world,0,0,0.0);
						}
						else
						{
							dWorldSetGravity (world,0,0,-9.8);
						}
						this->CenterOfMass(startpos);
						

					//dVector3 p = {40.0, 40.0, 2.0};
					//boxID = add_box(0.5, 80.0, 1.0, 2.0,p,false); 
					//dVector3 p2 = {40.0, -10.0, 2.0};
					//boxID = add_box(0.5, 80.0, 1.0, 2.0,p2,false);

					//dVector3 p3 = {0.0, 10.0, 2.0};
					//boxID = add_box(0.5, 1.0, 50.0, 2.0,p3,false); 
					//dVector3 p4 = {60.0, 10.0, 2.0};
					//boxID = add_box(0.5, 1.0, 50.0, 2.0,p4,false); 
					//	cout<<startpos[2]<<"\n";

						//joint2 = add_hinge(id2,id3,anchor2,yAxis,-0.8,0.8, 15.0); //was -1.4, 0.8
					}
				}
				else
				{
					dBodyID id1 = createShape(0, -1.0, NULL, finished);finished=false;//, false);  finished=false;
					dBodyID id2 = createShape(1, 0.0, NULL, finished);finished=false;//,false); 
					dBodyID id3 = createShape(2, 1.0, NULL, finished);//, false);
					dVector3 yAxis={0.0,-1.0,0.0};
					dVector3 anchor1 = {0,0,30};
					jointIDs[0] = add_hinge(id1,id2,anchor1,yAxis,-0.8,0.8, 15.0); //was -1.4, 0.8
					dVector3 anchor2 = {0,0,30};
					jointIDs[1] = add_hinge(id2,id3,anchor2,yAxis,-0.8,0.8, 15.0); //was -1.4, 0.8
	
					//add_fixed(id1, id2);
					//add_fixed(id2, id3);
					//connect first and last part with joints. do we need to make bodies of these parts to do that?
				}

			}
			//for (int i=0; i<NEAT::NumMachines; i++)
			//{
			//	createShape(xp);//(i+1)*3.0);
			//	partNumber = 0;

			//}
			//partNumber=NEAT::NumParts;
		}
		if (SHAPE_MODE==2)
		{
			float s=-1.0f;
			//cout<<"TIME: "<<timer<<"\n";
		
			//if (timer==250)
			//{
			//	this->CenterOfMass(startpos);
				//startpos[2]=1.0;
			//}

			//if (false)
			if (MOTORS_ENABLED)
			{
				for (int m=0; m<NEAT::NumMachines-1; m++)
				{
					if (m==0||m==NEAT::NumMachines-2)
					{
						if (m%2==0)
							actuateJoint(jointIDs[m], sin((float)tmpTimer));
							//dJointSetHingeParam(joints[ jointIDs[m] ],dParamVel,sin((float)tmpTimer)*5.0); //left knee	
						else
							actuateJoint(jointIDs[m], sin((float)tmpTimer-1.57f));
							//dJointSetHingeParam(joints[ jointIDs[m] ],dParamVel,sin((float)tmpTimer-1.57f)*5.0); //left knee	
					}
					else if (tmpTimer<1.0)
					{
						actuateJoint(jointIDs[m], sin((float)tmpTimer));
						//dJointSetHingeParam(joints[ jointIDs[m] ],dParamVel,sin((float)tmpTimer)*5.0); 
					}
					//s*=-1.0f;
				}
			}

			//Calcuate "Energy" consumption
			int c=0;
			//this->fitness+
			//vector<dVector3*> newPositions;
			double dist=0;
			//energyConsumption
			for (int m=0; m<NEAT::NumMachines; m++)
			{
				for (int x=0; x<machines[m].size(); x++)
				{
					const dReal *pos = dGeomGetPosition(machines[m][x]); //NEAT::NumParts/3.0
					//if (c==0)
				//	{
					//	cout<<pos[0]<<"\t"<<oldPosition[c][0]<<"\n";
				//	}
					dist = sqrt ( pow(pos[0]-(*oldPosition[c])[0],2) + pow(pos[1]-(*oldPosition[c])[1],2) +  pow(pos[2]-(*oldPosition[c])[2],2) );
					
					//cout<<dist<<"\n";
					energyConsumption += dist;
					
					//dVector3 newPos={pos[0],pos[1],pos[2]};
					//newPositions.push_back(&newPos);

					(*oldPosition[c])[0] = pos[0];//,pos[1],pos[2]};
					(*oldPosition[c])[1] = pos[1];//,pos[1],pos[2]};
					(*oldPosition[c])[2] = pos[2];//,pos[1],pos[2]};
					c++;
				}
			}
			//dist/=(NEAT::NumMachines*machines[0].size());
			//cout<<dist<<"\n";
			//energyConsumption/=(NEAT::NumMachines*machines[0].size());
			//Overwrite old positions
			//oldPosition.clear();
			//for (int i=0; i<newPositions.size(); i++)
			//	oldPosition.push_back(&newPositions[i]);

			
			if (WATER && false) 
			{
				float viscosity = 0.3;
				float AreaX=1.0;
				float AreaY =1.0f; 
				float AreaZ = 1.0f;
				for (int i=0; i<bodies.size(); i++)
				{
					const dReal *lvel = dBodyGetLinearVel(bodies[i]);
					const dReal *avel = dBodyGetAngularVel(bodies[i]);
					const dReal *R = dBodyGetRotation(bodies[i]);
					dReal ss[3];
					//dGeomBoxGetLengths(this->_id, ss);

					//dReal AreaX = ss[1] * ss[2];
					////dReal AreaY = ss[0] * ss[2];
				//	dReal AreaZ = ss[0] * ss[1];

					dReal nx = (R[0] * lvel[0] + R[4] * lvel[1] + R[8] * lvel[2]) *  AreaX;
					dReal ny = (R[1] * lvel[0] + R[5] * lvel[1] + R[9] * lvel[2]) * AreaY;
					dReal nz = (R[2] * lvel[0] + R[6] * lvel[1] + R[10] * lvel[2]) * AreaZ;

					dReal temp = -nx * viscosity;
					//cout<<"TEMP "<<temp<<"\n";

					dBodyAddForce(bodies[i], temp * R[0], temp * R[4], temp * R[8]);

					temp = -ny * viscosity;
					dBodyAddForce(bodies[i], temp * R[1], temp * R[5], temp * R[9]);

					temp =-nz * viscosity;
					dBodyAddForce(bodies[i], temp * R[2], temp * R[6], temp * R[10]);

					nx = (R[0] * avel[0] + R[4] * avel[1] + R[8] * avel[2]) * AreaZ;
					ny = (R[1] * avel[0] + R[5] * avel[1] + R[9] * avel[2]) * AreaX;
					nz = (R[2] * avel[0] + R[6] * avel[1] + R[10] * avel[2]) * AreaY;

					if (false)
					{
					temp = -nx * viscosity * 2;
					dBodyAddTorque(bodies[i], temp * R[0], temp * R[4], temp * R[8]);

					temp = -ny * viscosity * 2;
					dBodyAddTorque(bodies[i], temp * R[1], temp * R[5], temp * R[9]);

					temp = -nz * viscosity * 2;
					dBodyAddTorque(bodies[i], temp * R[2], temp * R[6], temp * R[10]);
					}
					//temp = -ny * viscosity;
					//dBodyAddForce(bodies[i], temp * R[1], temp * R[5], temp * R[9]);
				}
			}
			//!dJointSetHingeParam(joints[joint2],dParamVel,sin(1-(float)tmpTimer)*5.0); //left knee	
		}

		if (false && SHAPE_MODE==2)
		{
			const dReal *pos = dGeomGetPosition(machines[0][NEAT::NumParts/2.0]);
					
			dVector3 anchor1 = {pos[0],pos[1],pos[2]};

			const dReal *pos2 = dGeomGetPosition(machines[0][NEAT::NumParts/2.0-5]);
			dVector3 anchor2 = {pos2[0],pos2[1],pos2[2]};
						
			//{,0,30}; //{0,0,30}
			//joint1 = add_hinge(id1,id2,anchor1,yAxis,-0.8,0.8, 15.0); //was -1.4, 0.8
									
			////dMatrix3 RI;
			//dRSetIdentity (RI);
			//const dReal ss[3] = {1.02,1.02,1.02};

		//	startpos = dGeomGetPosition(machines[0][machines[0].size()/2]);

			//dWorldSetGravity (world,0,0,-9.8);
						//dsDrawBox (anchor1, RI, ss);
					//	dsDrawBox (anchor2, RI, ss);
		}

	//	int n1=NEAT::NumMachines/2;
		//int n2=machines[0].size()/2;


		
#ifdef GRAPHICS
								//dMass m2;f
		dMatrix3 RI;
		dRSetIdentity (RI);
		//	dVector3 v = {0.0+m.c[0], 0.0+m.c[1], 30.0+m.c[2]};
		const dReal ss2[3] = {0.5,0.5,0.5};
		//dVector3 pos={startx, starty, 0};
		//dsDrawBox (startpos, RI, ss2);
				dsSetColor(0.0, 0.0, 0.0);
		
		//dVector3 centermass;
		//CenterOfMassBody(centermass, tempbody[0]);
//		dMatrix3 RI;
		//dRSetIdentity (RI);
		//const dReal ss3[3] = {1.00,1.00,1.00};

		//dsDrawBox (centermass, RI, ss3);


		if (SHAPE_MODE==2)
		{
			for (int e=0; e<joints.size(); e++)
			{
				dVector3 pos;
				dJointGetHingeAnchor(joints[e], pos);
				//if (e==1)
				//{
					
				//	jointPosition.push_back(&pos);
				//}

				//const dReal *pos = dGeomGetPosition(machines[0][e]);
		
				//					dVector3 lengthVector = {0.0,0.0,0.5/2.0}; //{0.5,0.5,0.5};

				//dsDrawLine(jointPosition[e-1], jointPosition[e]);
		
							dVector3 newpos;

						//	dMULTIPLY0_331(newpos, dGeomGetOffsetRotation(machines[0][e]), lengthVector); //To world coordinates
							//dVector3 pos2 = {newpos[0]+pos[0],newpos[1]+pos[1],newpos[2]+pos[2]};
				dsDrawSphere(pos, RI, 0.3);//ss2);
			}
		//	for (int e=1; e<jointPosition.size(); e++)
		//	{
		//		dsDrawLine((*jointPosition[e-1]), (*jointPosition[e]));
		//	}
		}

		//cout<<startpos[2]<<"\n";
		if (timer%500==0)
			cout<<timer<<"\t"<<fitness()<<"\t"<<energyConsumption<<"\n";
#endif
		//cout<<timer<<"\n";
		tmpTimer+=0.03; //0.05  was 0.03
		timer++; 

	}

};



static dWorldID world;
static dSpaceID space;


static vector<dBodyID> bodies;
static vector<dGeomID> geoms;
static vector<Creature*> creatures;


static dJointGroupID contactgroup;

static bool show_contacts = true;

#define CYLRADIUS    0.6
#define CYLLENGTH    2.0
#define SPHERERADIUS 0.5


#ifdef dDOUBLE
#define dsDrawCapsule dsDrawCapsuleD
#define dsDrawCylinder dsDrawCylinderD
#define dsDrawSphere dsDrawSphereD
#define dsDrawBox dsDrawBoxD
#define dsDrawLine dsDrawLineD
#endif



// this is called by dSpaceCollide when two objects in space are
// potentially colliding.


static void nearCallback (void *data, dGeomID o1, dGeomID o2)
{
//	return;
	dBodyID test;
	//assert(o1);
	//assert(o2);
	
	dBodyID b1 = dGeomGetBody(o1);
	dBodyID b2 = dGeomGetBody(o2);

	
	//if(b1==b2) ODE is doing that
	//{
	//	cout<<"!!!!";
	//}

//	if (b1 && b2 && dAreConnectedExcluding (b1,b2,dJointTypeContact)) 
//	{
	//		return;
	//}

//	if (b1 && b2 && dAreConnected (b1,b2)) return;

	if(o1 == floorplane || o2 == floorplane)
	{
		if(o1==floorplane)
			test=b2;
		if(o2==floorplane)
			test=b1;
		//test should equal the body that is colliding with floor
		for(int x=0;x<creatures.size();x++)
		{
			int bsize=creatures[x]->bodies.size();
		//	for(int y=0;y<bsize;y++)
		//		if (test==creatures[x]->bodies[y])
			//		creatures[x]->onground[y]=true;
		}

	}
	else
		return;


	//else
	//{
	//	return; //1.24
	//}

//	else 
	//	return;

//	return;
	//Check if they are part of the same body	!
	const int N = 32;
	dContact contact[N];
	int n = dCollide (o1,o2,N,&(contact[0].geom),sizeof(dContact));
	//cout<<n<<"\n";
	if (n > 0) 
	{
				//dMatrix3 RI;
			//dRSetIdentity (RI);
			//const dReal ss[3] = {0.5,0.5,0.5};
		for (int i=0; i<n; i++) 
		{
			//if (o1 != floorplane && o2 != floorplane)
			//	dsDrawBox (contact[i].geom.pos,RI,ss);
			
			contact[i].surface.mode = 0;
			contact[i].surface.mu = dInfinity; //50.0; // was: dInfinity
			dJointID c = dJointCreateContact (world,contactgroup,&contact[i]);
			dJointAttach (c, dGeomGetBody(contact[i].geom.g1), dGeomGetBody(contact[i].geom.g2));
		}
	}
}


// start simulation - set viewpoint

static void start()
{
	srand(time(0));
	//-35.6, 9.9, 38.3  -29 -23.5 0
	#ifdef GRAPHICS
		static float xyz[3] = {-35.6,9.9,38.3};
		static float hpr[3] = {-29,-23.5,0};
		dsSetViewpoint (xyz,hpr);
	#endif
}

#ifdef GRAPHICS
// called when a key pressed
std::wstring PathCreator(const std::wstring& dir, wchar_t *fileName)
{
    return dir.substr( 0, dir.size() - 1 ) + fileName;
} 

std::vector< std::string > GetAllFiles()
{
    std::vector< std::string > allFlsArr;
    WIN32_FIND_DATA file;
	//SetFileApisToOem();
	char chFolderpath[_MAX_PATH];
  // CString::CString strExtension   = _T("*.txt");

	char *path=NULL;
	size_t size;
//	path=_getcwd(path,size);
	TCHAR pwd[MAX_PATH];
	GetCurrentDirectory(MAX_PATH,pwd);

	strcpy(chFolderpath, pwd);
   strcat(chFolderpath, _T("\\genomes\\"));
   //strcpy(chFolderpath, _T("c:\\genomes\\"));
   strcat(chFolderpath,  _T("*.txt"));

 //  cout<<chFolderpath<<"\n";
    HANDLE search_handle = FindFirstFile(chFolderpath, &file);
    if(search_handle)
    {
        do
        {
			// strcpy(chFolderpath, pwd);
			strcpy(chFolderpath, _T("genomes\\"));
			 //strcpy(chFolderpath, _T("c:\\genomes\\"));
		 strcat(chFolderpath,  file.cFileName);
		  
//		 cout<<chFolderpath<<"\n";
            //std::wstring path = PathCreator(dir,file.cFileName);
            //std::wcout << path << std::endl;
            allFlsArr.push_back(chFolderpath);
		//	cout<<"->\n";
		//	_tprintf (TEXT("  %s   \n"), file.cFileName);
			//cout << file.cFileName[0] << std::endl;
        }
        while(FindNextFile(search_handle, &file));
    }
  //  CloseHandle(search_handle);
	return allFlsArr;
}    

void loadFile(char file_name2[50])
{
	ifstream iFile(file_name2, ios::in);
	Genome* start_genome2 = new Genome(0, iFile);
	iFile.close();

	// Create CPPN from start genome
	cppn_network = start_genome2->genesis(0);

	char file_name5[50]; //overwriten this file will trigger a load
	sprintf(file_name5,"genome%d.txt",NEAT::REQUERY);
	cppn_network->genotype->print_to_filename(file_name5);
}


static void command (int cmd)
{
			std::vector<string> all_files;
	switch (cmd) 
	{
		case '1': NEAT::NumParts = 250;  RESET=true;

		break;
		case '2': 
			NEAT::NumParts = 500; RESET=true;
			break;
		case '3': 		NEAT::NumParts = 1000;  RESET=true; break;

		case '8':
			{
				static float xyz[3] = {-13.6, 3.3, 7.82};
			static float hpr[3] = {-28,-26.5,0};
			dsSetViewpoint(xyz, hpr);
				break;
			}
			//-13.979590 3.3818276 7.8299947 -28.500000 -26.000000 0
		
		case '9': 
			{
			static float xyz[3] = {-35.6,9.9,38.3};
			static float hpr[3] = {-29,-23.5,0};
			dsSetViewpoint(xyz, hpr);
				break;
			}
		case '0': //close
			{
			static float xyz[3]={-12, 6, 38} ;
			static float hpr[3] = {-22,-25,0};
			dsSetViewpoint(xyz, hpr);
				break;
			}
		case 'd': 
		//	NEAT::NumMachines = 1;
			STEP_MODE=!STEP_MODE; RESET=true;
		break;
		case 'm': MOTORS_ENABLED=!MOTORS_ENABLED; break;//switch motors on and off    saveToSTL(); break;
		case ']': NEAT::NumMachines ++; RESET=true; break;
		case '[': if (NEAT::NumMachines >0) NEAT::NumMachines --; RESET=true;  break;

			case 't': SAVE_STATE=true; break;
			case 'g': dWorldSetGravity (world,0,0,-9.8);break;
			case 'h': dWorldSetGravity (world,0,0,0.0); break;
	case 'e': 
		{
//		NUM_PRINTERS  = 1;

		char tmp[50];
		sprintf(tmp,"undue%d.txt",NEAT::REQUERY);
		cppn_network->genotype->print_to_filename(tmp);

		for (int e=1; e<10; e++)
		{

			Genome *gen = cppn_network->genotype->mutate();
			char file_name[50];
			sprintf(file_name,"genome%d.txt",e);
			gen->print_to_filename(file_name);
		}
		break;
		}
	case 'u': //undue
		{
		char tmp2[50];
		sprintf(tmp2,"undue%d.txt",NEAT::REQUERY);
		
		loadFile(tmp2);
		break;
		}
	case 's':
		//Genome *gen = cppn_network->genotype->mutate();
		char file_name[50];
		 time_t seconds;

		 seconds = time (NULL);
		sprintf(file_name,"genomes/genome%ld.txt",seconds);//NEAT::REQUERY);
		cppn_network->genotype->print_to_filename(file_name);
		cout<<"File saved\n";
		break;
	case ',': 
		//vector< std::wstring > all_genomes = 
		{
		all_files = GetAllFiles();
		FILE_COUNTER--;
//
		if (FILE_COUNTER<0) FILE_COUNTER=all_files.size()-1;

		cout<<all_files[FILE_COUNTER]<<"\n";//GetAllFiles()[0]; 
		
		Genome* start_genome;

		ifstream iFile(all_files[FILE_COUNTER], ios::in);
		start_genome = new Genome(0, iFile);
		iFile.close();

		// Create CPPN from start genome
		cppn_network = start_genome->genesis(0);

		char file_name4[50]; //overwriten this file will trigger a load
		sprintf(file_name4,"genome%d.txt",NEAT::REQUERY);
		cppn_network->genotype->print_to_filename(file_name4);

		break;
		}
	case '.': 
		{
		all_files = GetAllFiles();
		FILE_COUNTER++;
//
		if (FILE_COUNTER>all_files.size()-1) FILE_COUNTER=0;

		cout<<all_files[FILE_COUNTER]<<"\n";//GetAllFiles()[0]; 
		
		Genome* start_genome;

		ifstream iFile(all_files[FILE_COUNTER], ios::in);
		start_genome = new Genome(0, iFile);
		iFile.close();

		// Create CPPN from start genome
		cppn_network = start_genome->genesis(0);

		char file_name2[50]; //overwriten this file will trigger a load
		sprintf(file_name2,"genome%d.txt",NEAT::REQUERY);
		cppn_network->genotype->print_to_filename(file_name2);
		break;
		}
	case 'r':
		{
			all_files = GetAllFiles();
			remove(all_files[FILE_COUNTER].c_str());
			break;
		}
	case 'n': //create random new genome
		//NEAT::NumMachines  = 1;
		for (int e=1; e<10; e++) //was 10
		{
			Genome * start_genome = new Genome(0, CPPN_INPUTS, CPPN_OUTPUTS, CPPN_HIDDEN, CPPN_INPUTS+CPPN_OUTPUTS+CPPN_HIDDEN, false, LINK_PROB  );
			//new Genome(0, 1, 6, 2, 10, true, 0.7  );//Genome(id, iFile);
			// Create CPPN from start genome
			cppn_network = start_genome->genesis(0);

			char file_name3[50]; //overwriten this file will trigger a load
			sprintf(file_name3,"genome%d.txt",e);
			cppn_network->genotype->print_to_filename(file_name3);
		}
		break;
	case 'p': if (SHAPE_MODE==1) SHAPE_MODE=2; else SHAPE_MODE=1; break;
	case 'l': 
		char file_name2[50];
		sprintf(file_name2,"SAVED_genome%d.txt",NEAT::REQUERY);
		int id=0;
		Genome* start_genome2;

		cout<<"Loading file "<<file_name2<<"\n";
		ifstream iFile(file_name2, ios::in);
		start_genome2 = new Genome(id, iFile);
		iFile.close();

		// Create CPPN from start genome
		cppn_network = start_genome2->genesis(id);

		char file_name5[50]; //overwriten this file will trigger a load
		sprintf(file_name5,"genome%d.txt",NEAT::REQUERY);
		cppn_network->genotype->print_to_filename(file_name5);

		break;
	}
}

//_DEBUG;dSINGLE;WIN32;_CRT_SECURE_NO_DEPRECATE;DS_LIB;%(PreprocessorDefinitions)
void drawGeom (dGeomID g, const dReal *pos, const dReal *R, int show_aabb)
{
	int i;

	if (!g) return;
	if (!pos) pos = dGeomGetPosition (g);
	if (!R) R = dGeomGetRotation (g);



	int type = dGeomGetClass (g);
	if (type == dBoxClass) {
		dVector3 sides;
		dGeomBoxGetLengths (g,sides);
		dsDrawBox (pos,R,sides);
	}
	 else if (type == dConvexClass) 
    {
//#if 0
		//points2,
      dsDrawConvex(pos,R,planes,
		   planecount,
		   (float*)dGeomGetData(g),
		   pointcount,
		   polygons);
//#else
//      dsDrawConvex(pos,R,
//       Sphere_planes,
//		   Sphere_planecount,
//		   Sphere_points,
//		   Sphere_pointcount,
//		   Sphere_polygons);
//#endif
    }
	else if (type == dSphereClass) {
		dsDrawSphere (pos,R,dGeomSphereGetRadius (g));
	}
	else if (type == dCapsuleClass) {
		dReal radius,length;
		dGeomCapsuleGetParams (g,&radius,&length);
		dsDrawCapsule (pos,R,length,radius);
	}
	else if (type == dCylinderClass) {
		dReal radius,length;
		dGeomCylinderGetParams (g,&radius,&length);
		dsDrawCylinder (pos,R,length,radius);
	}
	else if (type == dGeomTransformClass) {
		dGeomID g2 = dGeomTransformGetGeom (g);
		const dReal *pos2 = dGeomGetPosition (g2);
		const dReal *R2 = dGeomGetRotation (g2);
		dVector3 actual_pos;
		dMatrix3 actual_R;
		dMULTIPLY0_331 (actual_pos,R,pos2);
		actual_pos[0] += pos[0];
		actual_pos[1] += pos[1];
		actual_pos[2] += pos[2];
		dMULTIPLY0_333 (actual_R,R,R2);
		drawGeom (g2,actual_pos,actual_R,0);
	}

//	show_aabb=true;

	if (show_aabb) {
		// draw the bounding box for this geom
		dReal aabb[6];
		dGeomGetAABB (g,aabb);
		dVector3 bbpos;
		for (i=0; i<3; i++) bbpos[i] = 0.5*(aabb[i*2] + aabb[i*2+1]);
		dVector3 bbsides;
		for (i=0; i<3; i++) bbsides[i] = aabb[i*2+1] - aabb[i*2];
		dMatrix3 RI;
		dRSetIdentity (RI);
		dsSetColorAlpha (1,0,0,0.5);
		dsDrawBox (bbpos,RI,bbsides);
	}
}
// render thingee

void render_body(dBodyID body, dGeomID geom)
{
	const dReal *CPos = dBodyGetPosition(body);
	const dReal *CRot = dBodyGetRotation(body);
	dReal cpos[3] = {CPos[0], CPos[1], CPos[2]};
	dReal * crot = new dReal[NUMBER_SEGMENTS*6];
	for (int i=0; i<NUMBER_SEGMENTS*6; i++)
		crot[i] = CRot[i];
	//= { CRot[0], CRot[1], CRot[2], CRot[3], CRot[4], CRot[5], CRot[6], CRot[7], CRot[8], CRot[9], CRot[10], CRot[11]}; //TODO 12

	if(dGeomGetClass(geom)==dCylinderClass)
	{
		dReal rad,length;
		dGeomCylinderGetParams(geom,&rad,&length);
		dsDrawCylinder
			(
			cpos,
			crot,
			length,
			rad
			); // single precision
	}
	else if(dGeomGetClass(geom)==dSphereClass)
	{
		dReal rad=dGeomSphereGetRadius(geom);

		dsDrawSphere
			(
			cpos,
			crot,
			rad
			); // single precision
	}
	else if(dGeomGetClass(geom)==dBoxClass)
	{
		dVector3 sides;
		dGeomBoxGetLengths(geom,sides);
		dsDrawBox
			(
			cpos,
			crot,
			sides
			);
	}
	delete[] crot;
}
#endif
// simulation loop

void simulationStep()
{
	double timestep=0.01; //0.01
	if(!bMoviePlay)
	{
		dVector3 v;
		dWorldGetGravity(world, v);
		if (v[2]!=0.0)
			dSpaceCollide (space,0,&nearCallback); //dWorldStep
		
		dWorldQuickStep(world, timestep);
	//	dWorldStep(world,timestep);//,10);
	}
	for(int x=0;x<creatures.size();x++)
		creatures[x]->Update(timestep);

	if(!bMoviePlay)
		dJointGroupEmpty (contactgroup);
}

static void create_world(Controller* controller,bool log=false)
{
	// create world
	//dRandSetSeed(dRand());
	dInitODE();


	world = dWorldCreate();
	space = dHashSpaceCreate (0);
	contactgroup = dJointGroupCreate (0);
	//was -9.8);//
	dWorldSetGravity (world,0,0,0.0); //-1.5, -9.8, -15.0
	floorplane = dCreatePlane (space,0,0,1, 0.0);
	dWorldSetERP(world,0.1); //0.1
	dWorldSetCFM(world,1E-4);

	Biped* biped = new Biped(log,bMoviePlay);
	dVector3 pos={0.0,0.0,0.0};
	biped->Create(world,space,pos,controller);

	//	Quadruped* quadruped = new Quadruped(log,bMoviePlay);
	//	dVector3  pos2={2.0,0.0,0.0};
	//	quadruped->Create(world,space,pos2,controller);

	//cout << "Total mass:" << biped->TotalMass() << endl;
	creatures.push_back(biped);

		
	srand(NEAT::REQUERY+ (unsigned)time( NULL ) );
	//	creatures.push_back(quadruped);
}

static void destroy_world()
{
			
	if (CPPN_DET_NUM_PARTS) //reset
		NEAT::NumParts=-1;

	dJointGroupEmpty (contactgroup);
	dJointGroupDestroy (contactgroup);

	for(int x=0;x<geoms.size();x++)
		dGeomDestroy(geoms[x]);

	for(int x=0;x<creatures.size();x++)
	{
		creatures[x]->Destroy();
		delete creatures[x];
	}
	creatures.clear();
	bodies.clear();
	geoms.clear();
	

	dSpaceDestroy (space);
	dWorldDestroy (world);
	dCloseODE();
}

#ifdef GRAPHICS
static void simLoop (int pause)
{
	static int timestep=0;

	if(bDisplay && (creatures[0]->abort()||timestep>1000) && bDoEvolution)
	{
		cout <<"Teriminating walker..." << endl;
		bDisplay=false;
		timestep=0;
		destroy_world();
		delete controller;
	}

	if(!bDisplay)
	{

		if(1)
		{
			biped_epoch(neatpop,Novelty);

			NEAT::Organism* champ;
			double max_fitness=0;
			std::vector<NEAT::Organism*>::iterator curorg;
			champ=*(neatpop->organisms.begin()); //Make sure at least something is chosen
			//Find the population champ
			for(curorg = neatpop->organisms.begin(); curorg != neatpop->organisms.end(); ++curorg) {
				if (((*curorg)->fitness>max_fitness)) {
					champ=(*curorg);
					max_fitness=champ->fitness;
				}
			}
			// Generate substrate genome from input file
			//char curword[100];
			//int id;
			//ifstream substrate_file(substrate_fn, ios::in);
			//substrate_file >> curword >> id;
			//Genome* substrate_genome = new Genome(id, substrate_file);

			// Create controller using generated substrate from champ CPPN
			//	controller = new CTRNNController(substrate_genome->genesis(id, champ->net, FIXED_LEG_SCALE));
			//controller = new CTRNNController(champ->net);

			//delete substrate_genome;
		}
		create_world(controller,true);
		cout << "Walker created..." << endl;
		bDisplay=true;
	}

	if (!pause) //&& !creatures[0]->abort())
	{
		// if(bSlowDown)
		//	usleep(10000);
		simulationStep();
		timestep+=1;  

		//if(timestep%100==0)
		//	cout << creatures[0]->fitness() << endl;
	} 

	dsSetColorAlpha (0.3,1.0,1.0,1.0);
	dsSetTexture(DS_WOOD);

	//if (timestep>30)
	//{
	for(int x=0;x<geoms.size();x++)
	{
		drawGeom(geoms[x],0,0,false);
	}


	dsSetTexture (DS_NONE);
		dMatrix3 RI;
	dRSetIdentity (RI);

		//	dVector3 sides;
		//dGeomBoxGetLengths (g,sides);
		//dsDrawBox (pos,R,sides);
	dsSetColorAlpha (0.0,0.0,1.0,0.1);

	for(int x=0;x<creatures.size();x++)
	{
		//int count=-2;
		//cout<<creatures[x]->geoms.size()<<"\n";

		int l= ((float)((NEAT::NumParts-creatures[x]->geoms.size())*0.5));
		dVector3 pos2 = {0, -l/2-0.25, 30}; 
		dVector3 sides = {0.4, l, 0.1};
	//	dGeomBoxGetLengths (g,sides);
		dsSetColor(0.0, 0.0, 1.0);
		//dsDrawBox (pos2,RI,sides);

		int left =(NEAT::NumParts-creatures[x]->geoms.size());

		//for (int st=0; st<

		for (int i=creatures[x]->geoms.size(); i<NEAT::NumParts; i++)
		{
			if (i%30==0)
			{
//				const dReal *pos = dGeomGetPosition (g);
				dsSetColor(1.0, 0.0, 0.0);
				dVector3 pos3 = {0, (i-creatures[x]->geoms.size())*(-0.5), 30}; 
				//cout <<(i-creatures[x]->geoms.size())<<"\n";
			//!	dsDrawSphere(pos3, RI, 0.4);
			}
		}
	//	cout<<"----\n";

		//for (int i=0; i<(150-creatures[x]->geoms.size()); i++)
		////{
		//}

		for(int y=0;y<creatures[x]->geoms.size();y++)
		{
			//dsSetColorAlpha (1.0,1.0,1.0,1.0);

			int start=0;
			//if (SHAPE_MODE==1) start = 2; //paint printer in different color

			if (y>=start && (y-start)%30==0)
			{
//				const dReal *pos = dGeomGetPosition (g);

				dsSetColor(1.0, 0.0, 0.0);
		//!		dsDrawSphere(dGeomGetPosition(creatures[x]->geoms[y]), RI, 0.4);
			}

			if (y>=start && y<NEAT::NumParts +2)
			{
				dsSetColor(0.0, 0.0, 1.0);
			}
			else if (y>NEAT::NumParts+3 && y <(NEAT::NumParts+2)*2)
			{
				dsSetColor(1.0, 0.0, 0.0);
			}
			else if (y>(NEAT::NumParts+2)*2+1)// && y <(NEAT::NumParts+2)*3)
			{
				dsSetColor(0.0, 1.0, 0.0);
			}
			else
			{
			//	dsSetColor(0.75, 0.75, 0.75);
				dsSetColor(1.0, 1.0, 1.0);
				//dsSetColorAlpha (0.75,0.75,0.75,0.3);

			}
			//dsSetColor(0.0, 1.0, 0.0);
			//count++;
			dsSetColor(0.0, 0.0, 1.0); //dsSetColor(1.0, 0.0, 1.0);
			dsSetColorAlpha (0.0,0.0,1.0,0.5);
			drawGeom(creatures[x]->geoms[y],0,0,false);
		}
	}
	//}
}
#endif

//For novelty search
void update_stuff(vector<float> &k, Creature* c,bool good=true)
{
	dVector3& o_com= ((Biped*)c)->orig_com;
	//dVector3& c_com= ((Biped*)c)->curr_com;
	dVector3 com;
	dVector3 delta;
	c->CenterOfMass(com);
	//calculate_delta(o_com,c_com,delta);
	calculate_delta(o_com,com,delta);
	//if(novelty_function==NF_COGSAMPSQ || 
	//	novelty_function == NF_COGSQ)
		calculate_power(delta,2);

	//if(novelty_function==NF_COGCUSAMP ||
	//	novelty_function == NF_COGCU)
	//	calculate_power(delta,3);

	//if(good)
	//{
		//cout<<delta[0]<<"\t"<<delta[1]<<"\n";
	k.push_back(delta[0]);
	k.push_back(delta[1]);
	//}
}

dReal evaluate_controller(Controller* controller,noveltyitem* ni=NULL,data_record* record=NULL,bool log=false)
{

	vector<float> k;
	dReal fitness;
	int timestep=0;
	const int simtime= 1500; //was 1500
	create_world(controller,log);
	while(!creatures[0]->abort() && timestep<simtime)
	{
		simulationStep();

		timestep+=1;

		if(timestep%100 == 0 && Novelty)//novelty_function % 2 == 1)
		{
			update_stuff(k,creatures[0]);
		}
	}

	//if(novelty_function%2==1)
	//{
	//	for(int x=timestep+1;x<=simtime;x++)
	//		if(x%100==0)
	//			update_stuff(k,creatures[0]); //,false);
	//}
	//else
	//{
	//	update_stuff(k,creatures[0]);
	//}

	//cout << timestep << endl;
	fitness=creatures[0]->fitness();
	//((Biped*)creatures[0])->lft.push_back(timestep);
	//((Biped*)creatures[0])->rft.push_back(timestep);
	
	if (ni!=NULL)
	{
		ni->novelty_scale = 1.0; //1.0 / avgDif;
		ni->data.push_back(k);
	}

	destroy_world();
	return fitness;
}


//display 1 NOOB02_F_120_4_2_bestgenome55_12.0736 4 120

//evolve novelty 4 120 output/run_150_10_2_ 2
//evolve output_ fitness 1 5 150
//symetry_seed.txt

Genome* start_genome;

int main(int argc, char **argv)
{
	// setup pointers to drawstuff callback functions(
	char curword[200];
	int id;

	first = true;

	//float f = 0.394893489348394893849238409284092384029384029348203949f;
	//cout<<f<<"\n";
	//return 0;
#ifdef GRAPHICS
	dsFunctions fn;
	fn.version = DS_VERSION;
	fn.start = &start;
	fn.step = &simLoop;

	fn.command = &command;
	fn.stop = 0;
	fn.path_to_textures = DRAWSTUFF_TEXTURE_PATH;
#endif 

	srand( (unsigned)time( NULL ) );
	load_neat_params("params.ne", true);
	dRandSetSeed(25);

	//load_parameters();
//	Genome* start_genome;

	//novelty_function seed hyperneat_variant
	//evolve output fitness bipedstartgenes substrate2i2df 7 0 1
	//display temp.dat substrate2i2df 0
	if(argc > 1)
	{
		if(strcmp(argv[1], "display") == 0)
		{

			NEAT::HYPERNEAT_SETTING = 0;
			NEAT::REQUERY = atoi(argv[2]); //right now stores the unique ID of this interactive window
			//cout<< (unsigned)time( NULL )<<"\n";



			srand( (unsigned)time( NULL ) );
		//	for (int i=0; i<10; i++)
			//	cout<<rand() <<"\n";

			//srand(24);
			//dRandSetSeed(35);//(unsigned)time( NULL )+ NEAT::REQUERY );
			
//			NUM_PARTS = atoi(argv[4]);
			//ifstream iFile(argv[2], ios::in);
			//ifstream substrate_file(argv[3], ios::in);
			//strcpy(substrate_fn, argv[3]);

			//cout<<"Reading in the start genome"<<endl;
			//Read in the start Genome
			//iFile >> curword;
			////iFile >> id;
		//	cout << "Reading in Genome id " << 0 << endl;
			//(int new_id,int i, int o, int n,int nmax, bool r, double linkprob)
			//int num_hidden=0;

			if (argc>3)
			{
				NEAT::NumMachines = atol(argv[4]);
			}
			if (argc>4)
			{
				NEAT::NumParts = atol(argv[5]);
			}
			
			
			if (NEAT::NumParts==-1) //add extra cppn output that determines the number of parts
			{
				CPPN_OUTPUTS = 3;
				CPPN_DET_NUM_PARTS=true;
				SEGMENTS_BETWEEN_MOTORS = NEAT::NumMachines; 
			}

			
			if (argc>3)
			{
				ifstream iFile(argv[3], ios::in);
				//ifstream substrate_file(argv[3], ios::in);
				//strcpy(substrate_fn, argv[3]);

				cout<<"Reading in the start genome"<<endl;
				//Read in the start Genome
				iFile >> curword;
				iFile >> id;
				cout << "Reading in Genome id " << id << endl;
				start_genome = new Genome(id, iFile);
				iFile.close();
			}
			else
			{
				start_genome = new Genome(0, CPPN_INPUTS, CPPN_OUTPUTS, CPPN_HIDDEN, CPPN_INPUTS+CPPN_OUTPUTS+CPPN_HIDDEN, false, LINK_PROB );//Genome(id, iFile);
			}


			//iFile.close();

			//if (argc==3)
			//	return 0;

			// Create CPPN from start genome
			cppn_network = start_genome->genesis(0);

			bSlowDown = true;
			bDisplay = true;
			bDoEvolution = false;
		}
		else
		{
			
			if(strcmp(argv[2], "novelty") == 0)
				Novelty = true;
			else
				Novelty = false;

			
			//if (argc>5)
			//{
				NEAT::NumMachines = atol(argv[3]);
			//}

			//if (argc>6)
			//{
				NEAT::NumParts = atol(argv[4]);
			//}

			if (NEAT::NumParts==-1)  //add extra cppn output that determines the number of parts
			{
				CPPN_OUTPUTS = 3;
				CPPN_DET_NUM_PARTS=true;
				SEGMENTS_BETWEEN_MOTORS = NEAT::NumMachines; //use nummachines as segments between motors if the cppn determines the num of parts
			}

			strcpy(outdir, argv[5]);




			//if (argc >= 5)
			//{
			//	strcpy(startgenes_fn, argv[4]);
			//	cout << "using " << startgenes_fn << " for start genes" << endl;
			//}

			if (argc >= 6)
			{
			//	strcpy(substrate_fn, argv[5]);
			//	cout << "using " << substrate_fn << " for substrate genes" << endl;
			}

			//if (argc >= 7)
			//{
			//atoi(argv[6]);
			//}
			novelty_function = NF_COGSAMPSQ;  //! TODO RISI  NF_FITSQSAMP;//

			//if (argc >= 7)
			//{

			//		dRandSetSeed(atol(argv[4]));
			NEAT::REQUERY = atol(argv[6]);

			
			srand(atol(argv[6])+ (unsigned)time( NULL ) );
				//srand(atol(argv[4]));//+(unsigned)time( NULL ) );



			if (argc>7)
			{
				ifstream iFile(argv[7], ios::in);
				//ifstream substrate_file(argv[3], ios::in);
				//strcpy(substrate_fn, argv[3]);

				cout<<"Reading in the start genome"<<endl;
				//Read in the start Genome
				iFile >> curword;
				iFile >> id;
				cout << "Reading in Genome id " << id << endl;
				start_genome = new Genome(id, iFile);
				iFile.close();
			}

			//}

			//if (argc >= 8)
			//{
			//	NEAT::HYPERNEAT_SETTING = atoi(argv[7]);
			//	cout <<"HyperNEAT variant "<<NEAT::HYPERNEAT_SETTING<<"\n";

			//}
			//if (argc >= 9)
			//{
			//	NEAT::REQUERY = atoi(argv[8]);
			//	cout <<"Requery CPPN "<<NEAT::REQUERY<<"\n";

			//}
		}	
	}


	if (bDoEvolution)
	{
		COLLISION_COUNTER=0;
		neatpop=biped_realtime(Novelty);
	//	cout<<COLLISION_COUNTER<<"\n";
		if(bNoVis)
		{
			int requery = NEAT::REQUERY; //backup
			//	NEAT::REQUERY = 0; //! initially don't requery either way

			for(int k=0; k < 501; k++)   //was 1000
			{
				biped_epoch(neatpop, Novelty);

			}
		}
		return 0;
	}

	// run simulation
#ifdef GRAPHICS
	if(!bNoVis)

	{

		if(!bDoEvolution)
			create_world(controller,true);
		
	
		//50, 80, 
		dsSimulationLoop(argc,argv,352, 288,&fn); // 352,288,

		cin.ignore();
		destroy_world();

	}
#endif
	//Determine number of segments that the CPPN outputs

					ifstream iFile(argv[3], ios::in);
				//ifstream substrate_file(argv[3], ios::in);
				//strcpy(substrate_fn, argv[3]);

				cout<<"Reading in the start genome"<<endl;
				//Read in the start Genome
				iFile >> curword;
				iFile >> id;
				cout << "Reading in Genome id " << id << endl;
				start_genome = new Genome(id, iFile);
				iFile.close();

				cppn_network = start_genome->genesis(0);


	double inputs[CPPN_INPUTS];
	inputs[0]=0.0; //pos
	inputs[1]=0; //distance from center
	//inputs[2]=1;//bias

	cppn_network->flush();
	cppn_network->load_sensors(inputs);
	//cout<<cppn_network->nodecount()<<"\n";
	for(int i = 0; i < CPPN_QUERIES; i++) //cppn_network->nodecount(); i++) //was 15
		cppn_network->activate();

	//also change number of motors
	NEAT::NumParts = (int)( cppn_network->outputs[2]->activation*240.0+60.0);//-cppn_network->outputs[3]->activation;

	//If the number of parts is determined by the CPPN the give number of machines give the number for a 80 piece robot. Therefore if the robot is longer we have to scale it
	NEAT::NumMachines = (int)(SEGMENTS_BETWEEN_MOTORS *(NEAT::NumParts/60.0));

	//make sure numparts is divisible by nummachines
	NEAT::NumParts= ((int)(NEAT::NumParts/NEAT::NumMachines))*NEAT::NumMachines;

	char file[50];
	FILE * pFile;
	sprintf(file,"cppn_logfile.txt",outdir);
	pFile = fopen (file,"a");
	if (pFile!=NULL)
	{
		char info[50];
		sprintf(info,"%f %f\n",NEAT::NumMachines,NEAT::NumParts);
		fputs (info, pFile);
		fclose (pFile);
	}

	cout << NEAT::NumMachines <<"\t"<<NEAT::NumParts<<"\n";


	cout << bNoVis << endl;
	if(bNoVis && !bDoEvolution)
	{
		return 0;
	}
}


NEAT::Population *biped_realtime(bool novelty) {
	NEAT::Population *pop;
	//NEAT::Genome *start_genome;
	char curword[20];
	int id=0;

	ostringstream *fnamebuf;
	int gen;

	double highscore;

	//ifstream iFile(startgenes_fn,ios::in);

//	std::string str = startgenes_fn;

	cout<<"Reading in the start genome"<<endl;
	//Read in the start Genome
	//iFile>>curword;
	//iFile>>id;
	if (start_genome==NULL)
	{
		cout<<"Reading in Genome id "<<id<<endl;

			//	ifstream iFile(argv[2], ios::in);
			//ifstream substrate_file(argv[3], ios::in);
			//strcpy(substrate_fn, argv[3]);

			cout<<"Reading in the start genome"<<endl;
			//Read in the start Genome
			//iFile >> curword;
			//iFile >> id;
			cout << "Reading in Genome id " << id << endl;
			//(int new_id,int i, int o, int n,int nmax, bool r, double linkprob)
			
			start_genome = new Genome(0, CPPN_INPUTS, CPPN_OUTPUTS, CPPN_HIDDEN, CPPN_INPUTS+CPPN_OUTPUTS+CPPN_HIDDEN, false, LINK_PROB  );//Genome(id, iFile);
	}
			//iFile.close();

	cout<<"Start Genome: "<<start_genome<<endl;

	//Spawn the Population from starter gene
	cout<<"Spawning Population off Genome"<<endl;

	pop=new NEAT::Population(start_genome,NEAT::pop_size);


	//Alternative way to start off of randomly connected genomes
	//pop=new Population(pop_size,7,1,10,false,0.3);

	cout<<"Verifying Spawned Pop"<<endl;
	pop->verify();

	//Start the evolution loop using rtNEAT method calls 
	biped_realtime_loop(pop, novelty);

	//if (str.find("_gen")==-1)
	//{
		delete start_genome;
	//}
	return pop;
}

void biped_epoch(NEAT::Population *pop,bool novelty) {

	vector<NEAT::Organism*>::iterator curorg;
	vector<NEAT::Species*>::iterator curspecies;

	vector<NEAT::Species*>::iterator curspec; //used in printing out debug info                                                         

	vector<NEAT::Species*> sorted_species;  //Species sorted by max fit org in Species                                                  

	ostringstream file_names;

	static int epoch=0;
	int pause;
	bool win=false;

	static double champ_fitness=0;

	NEAT::Organism *champ;

	ostringstream file_name;

	//Real-time evolution variables                                                                                             
	int offspring_count;
	NEAT::Organism *new_org;

	//We try to keep the number of species constant at this number                                                    
	int num_species_target=NEAT::pop_size/15;

	//This is where we determine the frequency of compatibility threshold adjustment
	int compat_adjust_frequency = NEAT::pop_size/10;
	if (compat_adjust_frequency < 1)
		compat_adjust_frequency = 1;


	//Rank all the organisms from best to worst in each species
	pop->rank_within_species();                                                                            

	//Assign each species an average fitness 
	//This average must be kept up-to-date by rtNEAT in order to select species probabailistically for reproduction
	pop->estimate_all_averages();

	epoch++;
	//if(epoch>1) exit(1);
	if(epoch%50==0)
	{
		char file[50];

		if(epoch%100==0)
		{
			sprintf(file,"%sgen%d",outdir,epoch);
			pop->print_to_file_by_species(file);

			sprintf(file,"%sarchive.dat",outdir);
			archive.Serialize(file);
		}


		//sprintf(file,"%srecord.dat",outdir);
		//Record.serialize(file);

		sprintf(file,"%sfittest",outdir);
		archive.serialize_fittest(file);


	}
	char file[50];
	FILE * pFile;
	sprintf(file,"%slogfile.txt",outdir);
	pFile = fopen (file,"a");
	if (pFile!=NULL)
	{
		char info[50];
		sprintf(info,"%d %f 0 0 0 0 0 0 0 0 0 0\n",epoch,champ_fitness);
		fputs (info, pFile);
		fclose (pFile);
	}

	cout << "Generation: " << epoch << " " << champ_fitness << " \n";

	//Now create offspring one at a time, testing each offspring,                                                               
	// and replacing the worst with the new offspring if its better
	for (offspring_count=0;offspring_count<NEAT::pop_size;offspring_count++) { //*10000

		if (offspring_count % NEAT::pop_size == 0)
		{

			//cout << "Generation: " << offspring_count/NEAT::pop_size << " " << champ_fitness << " \n";
			//***
			//float h=-1;
			//for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) {
			//	if ( (*curorg)->fitness > h) h = (*curorg)->fitness;
			//}
			//cout << "Best in current population "<<h<<"\n";
			//****
		}

		//Every pop_size reproductions, adjust the compat_thresh to better match the num_species_targer
		//and reassign the population to new species                                              
		if (offspring_count % compat_adjust_frequency == 0) {

			if(novelty)
			{	
				//update fittest individual list		
				archive.update_fittest(pop);
				//refresh generation's novelty scores
				archive.evaluate_population(pop,true);
			}
			int num_species = pop->species.size();
			double compat_mod=0.1;  //Modify compat thresh to control speciation                                                     

			// This tinkers with the compatibility threshold 
			if (num_species < num_species_target) {
				NEAT::compat_threshold -= compat_mod;
			}
			else if (num_species > num_species_target)
				NEAT::compat_threshold += compat_mod;

			if (NEAT::compat_threshold < 0.3)
				NEAT::compat_threshold = 0.3;

			//!			cout<<"compat_thresh = "<<NEAT::compat_threshold<<endl;

			//Go through entire population, reassigning organisms to new species                                                  
			for (curorg = (pop->organisms).begin(); curorg != pop->organisms.end(); ++curorg) {
				pop->reassign_species(*curorg);
			}
		}

		/*
		//For printing only
		for(curspec=(pop->species).begin();curspec!=(pop->species).end();curspec++) {
		cout<<"Species "<<(*curspec)->id<<" size"<<(*curspec)->organisms.size()<<" average= "<<(*curspec)->average_est<<endl;
		}

		cout<<"Pop size: "<<pop->organisms.size()<<endl;
		*/

		//Here we call two rtNEAT calls: 
		//1) choose_parent_species() decides which species should produce the next offspring
		//2) reproduce_one(...) creates a single offspring fromt the chosen species
		new_org=(pop->choose_parent_species())->reproduce_one(offspring_count,pop,pop->species);

		//Now we evaluate the new individual
		//Note that in a true real-time simulation, evaluation would be happening to all individuals at all times.
		//That is, this call would not appear here in a true online simulation.
		//cout<<"Evaluating new baby: "<<endl;

		data_record* newrec=new data_record();
		newrec->indiv_number=indiv_counter;

		//	new_org->noveltypoint = biped_evaluate(new_org,newrec);

		//***
		float fitness=0.0f;
		float current_f = 0.0f;

		float count=0;

		for (int i=0; i<1; i++)
		{



			//			LEG_SCALE = LEG_SIZES[i];
			count++;

			//LEG_ANGEL = 0.7;

			//!DENSITY = 4.5-LEG_SCALE;
			//2.5 = 2;
			//3.25 = 1.25;
			//4.0 = 0.5
			if (new_org->noveltypoint==NULL)
			{
				new_org->noveltypoint=(biped_evaluate(new_org,newrec));
				fitness += new_org->noveltypoint->fitness;
				current_f = new_org->noveltypoint->fitness;
			}
			else
			{
				data_record* temprec=new data_record();
				noveltyitem *new_item = biped_evaluate( new_org,temprec );
				fitness += new_item->fitness;

				current_f = new_item->fitness;
				delete temprec;

				//new_org->noveltypoint->data[0].insert(new_org->noveltypoint->data[0].end(), new_item->data[0].begin(), new_item->data[0].end() );
				//NEW NOVELTY METRIC
				for (int e=0; e<new_org->noveltypoint->data[0].size(); e++)
				{
					new_org->noveltypoint->data[0][e] += new_item->data[0][e];

					//TODO Also change metric in biped_realtime_loop
					//std::cout<<new_org->noveltypoint->data[0][e]<<"\t";
					//if ( abs(new_org->noveltypoint->data[0][e] - new_item->data[0][e]) < 0.1)
					//{
					//	//If the two values are to different make the behavior less novel
					//	new_org->noveltypoint->data[0][e] = 0.0; //or make them negative. e.g. * -0.1
					//}
					//else //average
					//{
					//	new_org->noveltypoint->data[0][e] = ( new_org->noveltypoint->data[0][e] + new_item->data[0][e])/2.0;
					//}
				}

				delete new_item;
			}

			//!was (*curorg)
			//if ( new_org->substrate_network!=NULL) delete new_org->substrate_network;


			//if (current_f<=5.0)
			//	break;
			//else
			//	fitness += 10.0;

			//			std::cout<<"EVALUATING "<<LEG_SCALE<<" fitness"<<(*curorg)->noveltypoint->fitness<<"\n";

			//			std::cout<<i<<" size "<<(*curorg)->noveltypoint->data[0].size()<<"\n";
		}

		//for (int e = 0; e<new_org->noveltypoint->data[0].size(); e++)
		//{
		//	new_org->noveltypoint->data[0][e]/=count;
		//}

		//new_org->noveltypoint->fitness = (float)( fitness/count);//was  fitness/NUMBER_TESTING_SIZES
		//******

		if (new_org->noveltypoint->fitness>champ_fitness) 
		{
			champ_fitness = new_org->noveltypoint->fitness;


			file_name.str("");
			file_name<<outdir<<"bestgenome"<<epoch<<"_"<<champ_fitness;
			new_org->gnome->print_to_filename(file_name.str().c_str()); 
#ifdef DEBUG_OUTPUT		
			file_name.str("");
			file_name<<outdir<<"bestgenome"<<offspring_count<<".xhtml";
			new_org->substrate_network->print_to_SVG(file_name.str().c_str());
#endif

		}

		new_org->noveltypoint->indiv_number = indiv_counter;
		//calculate novelty of new individual
		if(novelty)
		{
			archive.evaluate_individual(new_org,pop);
			new_org->fitness*=new_org->noveltypoint->novelty_scale;
			//RISI add again archive.update_fittest(new_org);
			//RISI Looks like this happens now anyway no matter if novelty or not
		}	
		else
		{
			new_org->fitness=new_org->noveltypoint->fitness;
		}
		archive.update_fittest(new_org); //TODO check if the fittest now also gets saved for fitness

		//!		if (new_org->substrate_network!=NULL) delete new_org->substrate_network;

		//add record of new indivdual to storage
		indiv_counter++;
		Record.add_new(newrec);

		//update fittest list

		if (win) {
			cout<<"WINNER"<<endl;
			//pop->print_to_file_by_species("rt_winpop");
			break;
		}


		//Now we reestimate the baby's species' fitness
		new_org->species->estimate_average();

		//Remove the worst organism                                                                                               
		pop->remove_worst();

	}

	if(novelty)
	{
		archive.end_of_gen_steady(pop);
		//archive.add_randomly(pop);
		archive.evaluate_population(pop,false);
		cout << "ARCHIVE SIZE:" << 
			archive.get_set_size() << endl;
	}

}


void biped_realtime_loop(NEAT::Population *pop,bool novelty) {
	vector<NEAT::Organism*>::iterator curorg;
	vector<NEAT::Species*>::iterator curspecies;

	vector<NEAT::Species*>::iterator curspec; //used in printing out debug info                                                         

	vector<NEAT::Species*> sorted_species;  //Species sorted by max fit org in Species                                                  

	int pause;
	bool win=false;

	//ostringstream file_name;
	//file_name.str("");
//	file_name<<startgenes_fn<<"_MIN_FITNESS";

	//std::ofstream fitness_log(file_name.str().c_str());

	double champ_fitness=-1.0;
	NEAT::Organism *champ;

	//Real-time evolution variables                                                                                             
	int offspring_count;
	NEAT::Organism *new_org;

	//We try to keep the number of species constant at this number                                                    
	int num_species_target=NEAT::pop_size/15;

	//This is where we determine the frequency of compatibility threshold adjustment
	int compat_adjust_frequency = NEAT::pop_size/10;
	if (compat_adjust_frequency < 1)
		compat_adjust_frequency = 1;

	//Initially, we evaluate the whole population                                                                               
	//Evaluate each organism on a test                               
	int max=0;
	float average=0;
	int org_count=0;

	float max_average_distance = 0.0f;
	float current_distance=-100;

	for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) {

		//shouldn't happen                                                                                                        
		if (((*curorg)->gnome)==0) {
			cout<<"ERROR EMPTY GENOME!"<<endl;
			cin>>pause;
		}

		float fitness=-10.0f;
		float current_f=0.0f;
		float count =0;

		double c = 1/2.7;
		//TODO do this for all different leg sizes
		int solvecount=0;

		org_count++;

		count++;

		(*curorg)->noveltypoint=(biped_evaluate((*curorg), NULL));
		//double temp = evaluate_controller(controller,NULL,NULL,true);
		//fitness_log<<c<<"\t"<<(*curorg)->noveltypoint->fitness<<"\t"<<(1/c)<<"\n";

		if ((*curorg)->noveltypoint->fitness>current_distance)
			current_distance = (*curorg)->noveltypoint->fitness;

		//current_distance+= (*curorg)->noveltypoint->fitness;

		//		cout<<c<<"\t"<<(*curorg)->noveltypoint->fitness<<"\t"<<(1/c)<<"\n";
		//if (temp > distance_dynamic) distance_dynamic = temp;

	//	if ( (*curorg)->substrate_network!=NULL) delete (*curorg)->substrate_network;
		//}
		//cout<<best<<"\n";
		//current_distance+=best;

		//distance_dynamic = -1.0;
		//c-=0.01;
		//if ((*curorg)->noveltypoint->fitness>=6.0) solvecount++;



		//}
		//current_distance/=13.0;
		if (current_distance > max_average_distance)
			max_average_distance = current_distance;

		average +=solvecount;
//		fitness_log<<"SOLVED "<<solvecount<<" average: "<<current_distance<<"\n";

		//	cout<<"SOLVED "<<solvecount<<" average: "<<current_distance<<"\n";
		if (solvecount > max) 
			max = solvecount;
		//}

	}
	average/=org_count;
	//fitness_log<<"MAX: "<<max<<" Average: "<<average<<" max average distance "<<max_average_distance<<"\n";
	//fitness_log<<"-------------------------------";
	//fitness_log.close();
	//exit(0);

	//Get ready for real-time loop
	if(novelty)
	{
		//assign fitness scores based on novelty
		archive.evaluate_population(pop,true);
		//add to archive
		archive.evaluate_population(pop,false);
	}

}

noveltyitem* biped_evaluate(NEAT::Organism *org,data_record* data)
{
char curword[100];
	int id = 0;

	noveltyitem *new_item = new noveltyitem;

	// Generate substrate genome from input file
	//RISI
	//if (substrate_genome!=NULL)
	//delete substrate_genome;

	//Load it from disk the first time

	// Store substrate's genome as the new novelty item's genotype
	//new_item->genotype = substrate_genome;
	new_item->genotype = new Genome(*org->gnome);

	// Create CPPN from the Organism's genome
	cppn_network = org->gnome->genesis((*org->net).net_id);

	new_item->phenotype = cppn_network;

	// Create the controller with the substrate Network
	//CTRNNController* cont = new CTRNNController(substrate_network);

		
			//Destroy();

		//bodies.clear();
		//geoms.clear();

		//Destroy();
		//Create(_worldID,_spaceID,pos,0);
		//partNumber=0;

	new_item->fitness = evaluate_controller(NULL, new_item, data);

	//if (creatures.size()==0)
	//	create_world(controller,false);

	//creatures[0]->resetWorld();
	//for (int i=0; i<100; i++)
	//	creatures[0]->Update(0.0f);
	
	//new_item->fitness = creatures[0]->fitness();
	//evaluate_controller(cont, new_item, data);

//	cout<< new_item->fitness <<"\n";

	// Clean up objects that aren't needed any more
	//delete cont;
	//	delete substrate_network;		risi. does this get delete somewhere else?

	//RISI now we want to delete it because we only load it at the beginning
	//	delete substrate_genome;

//destroy_world();

	return new_item;

}
