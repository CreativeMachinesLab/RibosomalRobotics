#ifndef _QUADRUPED
#define _QUADRUPED

#include <vector>
#include <ode/odeconfig.h>
#include <iostream>
#include "Controller.h"
#include "ode/ode.h"

class Quadruped : public Creature
{
public:
//	static dReal D_CONSTANT= 0.0;

	//NEW
	dJointID joints[8];
	static dBodyID body[9];

	int step;
	dVector3 orig_com;
	dVector3 orig_left;
	dVector3 orig_right;
	dVector3 curr_com;
	bool log;
	std::ofstream* logfile;
	dJointFeedback feedback[6];



	//keeping track of foot positionxorz
	bool leftdown;
	bool rightdown;

	bool leftrigid;
	bool rightrigid;

	int lastdown; //which foot was last down

	std::vector<float> lft; //left foot time
	std::vector<float> lfx; //x
	std::vector<float> lfy; //y
	std::vector<float> rft; //right foot time
	std::vector<float> rfx; //x
	std::vector<float> rfy; //y

	Quadruped(bool logging=false,bool movie=false):Creature(logging,movie) {
		step=0;

		leftdown=false;
		rightdown=false;

		leftrigid=false;
		rightrigid=false;

		lastdown=0;

		log=logging;

		for(int x=0;x<6;x++)
		{
		//!	p_terms.push_back(P_CONSTANT);
		//!	d_terms.push_back(D_CONSTANT);

			//desired_angles.push_back(0.0);
			current_angles.push_back(0.0);
			delta_angles.push_back(0.0);
			desired_angvel.push_back(0.0);
			lo_limit.push_back(0.0);
			hi_limit.push_back(0.0);
		}

		for(int x=0;x<8;x++)
			sensors.push_back(0.0);
	} 
//	int add_foot(dReal density, dReal radius, const dVector3 p) {}
	virtual dReal fitness()
	{
		return 0.0;
	}

	~Quadruped()
	{
		//if(log)
		//	delete logfile;
	}
	void CreateBox( int index, double x, double y, double z, double length, double width, double height)
	{
	}

	void createHing(int index, dBodyID firstBody, dBodyID secondBody)
	{
		
	///	joints[index] = dJointCreateHinge(world,0);
	//	dJointAttach(joints[index],firstBody,secondBody);
	//	dJointSetHingeAnchor(joints[index],jointPosition[0],jointPosition[1],jointPosition[2]);
	//	dJointSetHingeAxis(joints[index],jointNormal[0],jointNormal[1],jointNormal[2]);
	}

	virtual bool abort() { return false;}
	void create_leg(dVector3 offset);
	virtual void Create(dWorldID worldi, dSpaceID spacei, dVector3 posi,Controller* cont);
//	void print_behavior() {}
	// 2 input
	virtual void Update(double timestep) {}
	// 6 input
	virtual void Update_dnu(double timestep) {}
};

#endif
