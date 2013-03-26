#ifndef _CONTROLLER
#define _CONTROLLER

#include <vector>
#include <ode/odeconfig.h>
#include <iostream>
#include <ode/odeconfig.h>
#include "ode/ode.h"
//#include <windows.h>

class Controller
{
public:
	std::vector<dReal> outs;
	int size;
	bool scale;
	bool debug;
	Controller() { } 
	Controller(int s,bool deb=false)
	{
		debug=deb;
		scale=true;
		size=s;
		for(int x=0;x<size;x++)
			outs.push_back(0.0);
	}
	virtual void update(double time,std::vector<dReal> sensors) {}
	virtual std::vector<dReal>* get_outputs()
	{
		return &outs;
	}
	virtual void reinforce(double r) {}
	virtual void revertChanges() {}
	virtual void resetCTRNN() {}
	virtual ~Controller() {}
};

class Creature
{
public:
//	static bool INItIAL_COLLESION;

	bool movie_rec;
	bool movie_play;

	Controller* controller;
	std::vector<dGeomID> geoms;
	std::vector<dBodyID> bodies;
	std::vector<bool> onground;
	std::vector<dJointID> joints;

	std::vector<dReal> current_angles;
//!	vector<dReal> desired_angles;

	std::vector< std::vector<dGeomID> > machines;

	std::vector<dReal> lo_limit;
	std::vector<dReal> hi_limit;

	std::vector<dReal> sensors;

	//vector<dReal> sensors;
	std::vector<dReal> desired_angvel;
	std::vector<dReal> delta_angles;

	std::vector<dReal> p_terms;
	std::vector<dReal> d_terms;
	
	int startx,starty;
	dVector3 startpos;
	dWorldID world;
	dSpaceID space;
	dVector3 pos;

	Creature(bool mov=false,bool play=false) {
		if(play)
			mov=false;
		movie_play=play;
		movie_rec=mov;
	}

	virtual void Update(double timestep) {}
	
	dReal TotalMass()
	{
		dReal total_mass=0.0;
		for(int x=0;x<bodies.size();x++)
		{
			dMass m;
			dBodyGetMass(bodies[x],&m);
			total_mass+=m.mass;
		}
		return total_mass;
	}

	void CenterOfMassBody(dVector3 center, dBodyID id)
	{
		dReal total_mass=0.0;
		dVector3 accum={0.0,0.0,0.0};
		const dReal* bpos;
	//	for(int x=0;x<bodies.size();x++)
	//	{
			dMass m;
			dBodyGetMass(id,&m);
			bpos=dBodyGetPosition(id);
			total_mass+=m.mass;
			for(int y=0;y<3;y++)
				accum[y]+=m.mass*bpos[y];
	//	}

		for(int x=0;x<3;x++)
			center[x]=accum[x]/total_mass;
	}

	void CenterOfMass(dVector3 center)
	{
		dReal total_mass=0.0;
		dVector3 accum={0.0,0.0,0.0};
		const dReal* bpos;
		for(int x=0;x<bodies.size();x++)
		{
			dMass m;
			dBodyGetMass(bodies[x],&m);
			bpos=dBodyGetPosition(bodies[x]);
			total_mass+=m.mass;
			for(int y=0;y<3;y++)
				accum[y]+=m.mass*bpos[y];
		}

		for(int x=0;x<3;x++)
			center[x]=accum[x]/total_mass;
	}

	virtual void Create(dWorldID worldi, dSpaceID spacei, dVector3 posi,Controller* cont)
	{
		controller=cont;
		world=worldi;
		space=spacei;
		pos[0]=posi[0];
		pos[1]=posi[1];
		pos[2]=posi[2];
	}

	virtual void Destroy()
	{
		for(int x=0;x<geoms.size();x++)
			dGeomDestroy(geoms[x]);
	}



	virtual bool abort()=0;
	virtual dReal fitness()=0;
	int add_fixed(int b1,int b2)
	{
		dBodyID bd1=bodies[b1];
		dBodyID bd2;
		if (b2!=(-1))
		{
			bd2=bodies[b2];
		}
		else
		{
			bd2=0;
		}
		dJointID tempjoint = dJointCreateFixed(world,0);
		dJointAttach(tempjoint,bd1,bd2);
		dJointSetFixed(tempjoint);
		joints.push_back(tempjoint);
		return joints.size()-1;
	}

	int add_fixed(dBodyID bd1,dBodyID bd2)
	{
		dJointID tempjoint = dJointCreateFixed(world,0);
		dJointAttach(tempjoint,bd1,bd2);
		dJointSetFixed(tempjoint);
		joints.push_back(tempjoint);
		return joints.size()-1;
	}

	int add_hinge(dBodyID b1,dBodyID b2,dVector3 anchor,dVector3 axis,dReal lostop, dReal histop, dReal fmax)
	{
		dJointID tempjoint = dJointCreateHinge(world,0);
		dJointAttach(tempjoint,b1,b2);
		dJointSetHingeAnchor(tempjoint,pos[0]+anchor[0],pos[1]+anchor[1],pos[2]+anchor[2]);
		dJointSetHingeAxis(tempjoint,axis[0],axis[1],axis[2]);
		dJointSetHingeParam(tempjoint,dParamLoStop,lostop);
		dJointSetHingeParam(tempjoint,dParamHiStop,histop);
		dJointSetHingeParam(tempjoint,dParamFMax,fmax);
		joints.push_back(tempjoint);

		return joints.size()-1;
	}

	int add_hinge(int b1,int b2,dVector3 anchor,dVector3 axis,dReal lostop, dReal histop, dReal fmax)
	{
		dJointID tempjoint = dJointCreateHinge(world,0);
		dJointAttach(tempjoint,bodies[b1],bodies[b2]);
		dJointSetHingeAnchor(tempjoint,pos[0]+anchor[0],pos[1]+anchor[1],pos[2]+anchor[2]);
		dJointSetHingeAxis(tempjoint,axis[0],axis[1],axis[2]);
		dJointSetHingeParam(tempjoint,dParamLoStop,lostop);
		dJointSetHingeParam(tempjoint,dParamHiStop,histop);
		dJointSetHingeParam(tempjoint,dParamFMax,fmax);
		joints.push_back(tempjoint);

		return joints.size()-1;
	}

	int add_universal(int b1,int b2,dVector3 anchor, dVector3 axis1, dVector3 axis2, dReal lostop1, dReal histop1, dReal lostop2, dReal histop2, dReal fmax1, dReal fmax2)
	{
		dJointID tempjoint = dJointCreateUniversal(world,0);
		dJointAttach(tempjoint,bodies[b1],bodies[b2]);
		dJointSetUniversalAnchor(tempjoint,pos[0]+anchor[0],pos[1]+anchor[1],pos[2]+anchor[2]);
		dJointSetUniversalAxis1(tempjoint,axis1[0],axis1[1],axis1[2]);
		dJointSetUniversalAxis2(tempjoint,axis2[0],axis2[1],axis2[2]);
		dJointSetUniversalParam(tempjoint,dParamLoStop,lostop1);
		dJointSetUniversalParam(tempjoint,dParamHiStop,histop1);
		dJointSetUniversalParam(tempjoint,dParamLoStop2,lostop2);
		dJointSetUniversalParam(tempjoint,dParamHiStop2,histop2);
		dJointSetUniversalParam(tempjoint,dParamFMax,fmax1);


		dJointSetUniversalParam(tempjoint,dParamFMax2,fmax2);

		joints.push_back(tempjoint);
		return joints.size()-1;
	}
	int add_box(dReal density, dReal lx, dReal ly, dReal lz, const dVector3 p,bool addToBodies=true)
	{
		dBodyID tempbody;
		dGeomID tempgeom;
		dMass m;
		tempbody = dBodyCreate (world);
		dMassSetBoxTotal(&m,density,lx,ly,lz);
		tempgeom = dCreateBox(0,lx,ly,lz);
		dGeomSetBody(tempgeom,tempbody);
		dBodySetPosition(tempbody,pos[0]+p[0],pos[1]+p[1],pos[2]+p[2]);
		dSpaceAdd(space,tempgeom);

		if (addToBodies)
		{
			bodies.push_back(tempbody);
			onground.push_back(false);
		}
		geoms.push_back(tempgeom);
		return bodies.size()-1;
	}

	int add_sphere(dReal density, dReal radius, const dVector3 p)
	{

		dBodyID tempbody;
		dGeomID tempgeom;
		dMass m;
		tempbody = dBodyCreate (world);
		dMassSetSphereTotal(&m,density,radius);
		tempgeom = dCreateSphere(0,radius);
		dGeomSetBody(tempgeom,tempbody);
		dBodySetMass(tempbody,&m);
		dBodySetPosition(tempbody,pos[0]+p[0],pos[1]+p[1],pos[2]+p[2]);
		dSpaceAdd(space,tempgeom);

		bodies.push_back(tempbody);
		geoms.push_back(tempgeom);
		onground.push_back(false);
		return bodies.size()-1;

	}

	int add_cylinder(dReal a[], dReal density,dReal length, dReal radius, const dVector3 p,dBodyID* k=NULL)
	{
		//dReal a[]={0.0,0.0,0.0};

		dBodyID tempbody;
		dGeomID tempgeom;
		dMass m;

		tempbody = dBodyCreate (world);
		if(k!=NULL)
			(*k)=tempbody;
		dQuaternion q;
		//if (axis==1)
		//{
		//	a[1]=1.0;
		//}
		//else if (axis==2)
		//{
		//	a[0]=1.0;
		//}
		//else
		//{
		//	a[2]=1.0;
		//}
		dMatrix3 R;
		//dRFromEulerAngles(R,a[0],a[1],a[2]);
		dRFromZAxis(R,a[0],a[1],a[2]);
		dBodySetRotation(tempbody,R);

	//	dQFromAxisAndAngle (q,a[0],a[1],a[2], M_PI * 0.5);
	//	dBodySetQuaternion (tempbody,q);
		dMassSetCylinderTotal (&m,density,3,radius,length);
		dBodySetMass (tempbody,&m);
		tempgeom = dCreateCylinder(0, radius, length);
		dGeomSetBody (tempgeom,tempbody);
		dBodySetPosition (tempbody, pos[0]+p[0],pos[1]+p[1], pos[2]+p[2]);		
		dSpaceAdd (space, tempgeom);

		geoms.push_back(tempgeom);
		bodies.push_back(tempbody);
		onground.push_back(false);
		return bodies.size()-1;
	}


	int add_cylinder(int axis, dReal density,dReal length, dReal radius, const dVector3 p,dBodyID* k=NULL)
	{
		dReal a[]={0.0,0.0,0.0};

		dBodyID tempbody;
		dGeomID tempgeom;
		dMass m;

		tempbody = dBodyCreate (world);
		if(k!=NULL)
			(*k)=tempbody;
		dQuaternion q;
		if (axis==1)
		{
			a[1]=1.0;
		}
		else if (axis==2)
		{
			a[0]=1.0;
		}
		else
		{
			a[2]=1.0;
		}
		dQFromAxisAndAngle (q,a[0],a[1],a[2], M_PI * 0.5);
		dBodySetQuaternion (tempbody,q);
		dMassSetCylinderTotal (&m,density,axis,radius,length);
		dBodySetMass (tempbody,&m);
		tempgeom = dCreateCylinder(0, radius, length);
		dGeomSetBody (tempgeom,tempbody);
		dBodySetPosition (tempbody, pos[0]+p[0],pos[1]+p[1], pos[2]+p[2]);		
		dSpaceAdd (space, tempgeom);

		geoms.push_back(tempgeom);
		bodies.push_back(tempbody);
		onground.push_back(false);
		return bodies.size()-1;
	}
	virtual ~Creature() {}
};

#endif
