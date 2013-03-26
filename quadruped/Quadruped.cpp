#include <ode/odeconfig.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include "neat.h"
#include "organism.h"
#include "datarec.h"
#include "Quadruped.h"

#define GRAPHICS

using namespace std;

static dReal SCALE_FACTOR = 2;
static dReal FOOTX_SZ =1.0/SCALE_FACTOR;
static dReal FOOTY_SZ =0.5/SCALE_FACTOR;
static dReal FOOTZ_SZ =1.0/SCALE_FACTOR;
static dReal LLEG_LEN =1.0/SCALE_FACTOR;
static dReal LLEG_RAD =0.2/SCALE_FACTOR;
static dReal ULEG_LEN =1.0/SCALE_FACTOR;
static dReal ULEG_RAD =0.2/SCALE_FACTOR;
static dReal TORSO_LEN =1.0/SCALE_FACTOR;
static dReal TORSO_RAD =0.3/SCALE_FACTOR;
static dReal ORIG_HEIGHT= (TORSO_RAD/2.0+ULEG_LEN+LLEG_LEN+FOOTZ_SZ);
static dReal DENSITY=0.5;
static dReal TORSO_DENSITY=1.0;
static dReal FOOT_DENSITY=0.1;
static dReal MAXTORQUE_FOOT= 10.0;
static dReal MAXTORQUE_KNEE= 5.0;
static dReal MAXTORQUE_HIPMINOR= 5.0;
static dReal MAXTORQUE_HIPMAJOR= 5.0;
static dReal P_CONSTANT= 9.0;
static dReal FOOTFACTOR= 5.0;


void Quadruped::create_leg(dVector3 offset)
{
	dVector3 xAxis={1.0,0.0,0.0};
	dVector3 yAxis={0.0,-1.0,0.0};
	dVector3 zAxis={0.0,0.0,1.0};

	dVector3 p={offset[0],offset[1],offset[2]};

	//dVector3 foot_pos = {p[0]-0.3*FOOTX_SZ,p[1],p[2]+(FOOTZ_SZ/2.0)};
	//int foot=add_box(DENSITY,FOOTX_SZ,FOOTY_SZ,FOOTZ_SZ,foot_pos);

	dVector3 foot_pos = {p[0],p[1],p[2]+(FOOTZ_SZ/2.0)};

	//int foot=add_sphere(FOOT_DENSITY,FOOTZ_SZ/2.0,foot_pos);


	//dVector3 lower_pos = {p[0],p[1],p[2]+FOOTZ_SZ+LLEG_LEN/2.0};
	//int lowerleg = add_cylinder(3,DENSITY,LLEG_LEN,LLEG_RAD,lower_pos);
	//dVector3 upper_pos = {p[0],p[1],p[2]+FOOTZ_SZ+LLEG_LEN+ULEG_LEN/2.0};
	//int upperleg = add_cylinder(3,DENSITY,ULEG_LEN,ULEG_RAD,upper_pos);

	double M_BODY_LEN = 0.7;

	dVector3 lower_pos = {p[0]+LLEG_LEN,p[1],p[2]};
	int lowerlegLeft = add_cylinder(1,DENSITY,LLEG_LEN,LLEG_RAD,lower_pos);
	
	dVector3 upper_pos = {p[0]+LLEG_LEN+ULEG_LEN,p[1],p[2]};
	int upperlegLeft = add_cylinder(1,DENSITY,ULEG_LEN,ULEG_RAD,upper_pos);
//----

	dVector3 lower_pos_right = {p[0]-LLEG_LEN,p[1],p[2]};			//	dVector3 lower_pos_right = {-M_BODY_LEN/2-p[0]-LLEG_LEN,p[1],p[2]};
	int lowerlegRight = add_cylinder(1,DENSITY,LLEG_LEN,LLEG_RAD,lower_pos_right);
	
	dVector3 upper_pos_right = {p[0]-LLEG_LEN-ULEG_LEN,p[1],p[2]};
	int upperlegRight = add_cylinder(1,DENSITY,ULEG_LEN,ULEG_RAD,upper_pos_right);
//-----

	dVector3 foot_joint1 = {M_BODY_LEN/2+p[0]+LLEG_LEN,p[1],p[2]};
	//dVector3 knee_joint_a = {p[0]+ULEG_LEN+LLEG_LEN,p[1],p[2]};
	
	dVector3 foot_joint2 = {-M_BODY_LEN/2-p[0]-LLEG_LEN,p[1],p[2]};
	//dVector3 knee_joint_a = {p[0]+ULEG_LEN+LLEG_LEN,p[1],p[2]};

	dVector3 body_joint1 = {p[0],p[1],p[2]};
	dVector3 body_joint2 = {p[0],p[1],p[2]};

	//CREATE THE MAIN BODY
	dVector3 body_p = {p[0],p[1],p[2]}; //-M_BODY_LEN/2	-M_BODY_LEN/2
	int mainBody = add_box(DENSITY, M_BODY_LEN, M_BODY_LEN, 0.1, body_p);

	//add_fixed(foot,lowerleg);
	//add_universal(foot,lowerleg,foot_joint_a,xAxis,yAxis,-0.1,0.1,-0.1,0.1,MAXTORQUE_FOOT,MAXTORQUE_FOOT);

	add_hinge(lowerlegLeft,upperlegLeft,foot_joint1,yAxis,-1.4,0.0,MAXTORQUE_KNEE);
	add_hinge(lowerlegRight,upperlegRight,foot_joint2,yAxis,-1.4,0.0,MAXTORQUE_KNEE);
	
	add_hinge(upperlegLeft,mainBody,body_joint1,yAxis,-1.4,0.0,MAXTORQUE_KNEE);
	add_hinge(upperlegRight,mainBody,body_joint2,yAxis,-1.4,0.0,MAXTORQUE_KNEE);
}


void Quadruped::Create(dWorldID worldi, dSpaceID spacei, dVector3 posi,Controller* cont)
{
	
	Creature::Create(worldi,spacei,posi,cont);

	double M_BODY_LEN = 0.7;
	
	dVector3 xAxis={1.0,0.0,0.0};
	dVector3 yAxis={0.0,-1.0,0.0};
	dVector3 zAxis={0.0,0.0,1.0};
//	double thickness

	//CREATE THE MAIN BODY
	dVector3 body_p = {0,0,LLEG_LEN}; //-M_BODY_LEN/2	-M_BODY_LEN/2
	int mainBody = add_box(DENSITY, M_BODY_LEN, M_BODY_LEN, 0.1, body_p);

	//CREATE LEG 3 	
	dVector3 upper_pos = {M_BODY_LEN/2+LLEG_LEN/2,0.0,0.5};
	int upperlegLeft = add_cylinder(1,DENSITY,ULEG_LEN,ULEG_RAD,upper_pos);

	dVector3 lower_pos = {M_BODY_LEN/2+LLEG_LEN,0.0,0.0};
	int lowerlegLeft = add_cylinder(3,DENSITY,LLEG_LEN,LLEG_RAD,lower_pos);


	dVector3 foot_joint1 = {M_BODY_LEN/2+LLEG_LEN,0.0,0.5};

	add_hinge(lowerlegLeft,upperlegLeft,foot_joint1,yAxis,-1.4,0.0,MAXTORQUE_KNEE);

	dVector3 body_joint1 = {M_BODY_LEN/2,0.0,0.5-0.1};
	add_hinge(upperlegLeft,mainBody,body_joint1,yAxis,-1.4,0.0,MAXTORQUE_KNEE);
	return;

//----
	//CREATE LEG 1

	dVector3 lower_pos_right = {-LLEG_LEN,0.0,0.0};			//	dVector3 lower_pos_right = {-M_BODY_LEN/2-p[0]-LLEG_LEN,p[1],p[2]};
	int lowerlegRight = add_cylinder(1,DENSITY,LLEG_LEN,LLEG_RAD,lower_pos_right);
	
	dVector3 upper_pos_right = {0.0-LLEG_LEN-ULEG_LEN,0.0,0.0};
	int upperlegRight = add_cylinder(1,DENSITY,ULEG_LEN,ULEG_RAD,upper_pos_right);
		
	dVector3 foot_joint2 = {-M_BODY_LEN/2-0.0-LLEG_LEN,0.0,0.0};
	add_hinge(lowerlegRight,upperlegRight,foot_joint2,yAxis,-1.4,0.0,MAXTORQUE_KNEE);
	
	dVector3 body_joint2 = {0.0,0.0,0.0};
	
	add_hinge(upperlegRight,mainBody,body_joint2,yAxis,-1.4,0.0,MAXTORQUE_KNEE);

//-----


	//dVector3 knee_joint_a = {p[0]+ULEG_LEN+LLEG_LEN,p[1],p[2]};

	//dVector3 knee_joint_a = {p[0]+ULEG_LEN+LLEG_LEN,p[1],p[2]};




	//add_fixed(foot,lowerleg);
	//add_universal(foot,lowerleg,foot_joint_a,xAxis,yAxis,-0.1,0.1,-0.1,0.1,MAXTORQUE_FOOT,MAXTORQUE_FOOT);

}

