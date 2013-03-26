#include "network.h"
#include <string.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>

using namespace NEAT;

void Network::ctrnn_dynamics()
{
	std::ofstream out("out.dat");
	init_ctrnn();
	for(int x=0;x<300;x++)
	{
		activate_ctrnn(0.02);
		for(int y=0;y<6;y++)
			out << outputs[y]->output << " "; 
		out << std::endl;
	}

}


Network::Network(std::vector<NNode*> in,std::vector<NNode*> out,std::vector<NNode*> all,int netid) {
	inputs=in;
	outputs=out;
	all_nodes=all;
	name=0;   //Defaults to no name  ..NOTE: TRYING TO PRINT AN EMPTY NAME CAN CAUSE A CRASH
	numnodes=-1;
	numlinks=-1;
	net_id=netid;
	adaptable=false;

	//fnetout = new std::ofstream("act.txt");

}

Network::Network(std::vector<NNode*> in,std::vector<NNode*> out,std::vector<NNode*> all,int netid, bool adaptval) {
	inputs=in;
	outputs=out;
	all_nodes=all;
	name=0;   //Defaults to no name  ..NOTE: TRYING TO PRINT AN EMPTY NAME CAN CAUSE A CRASH                                    
	numnodes=-1;
	numlinks=-1;
	net_id=netid;
	adaptable=adaptval;
}


Network::Network(int netid) {
	name=0; //Defaults to no name
	numnodes=-1;
	numlinks=-1;
	net_id=netid;
	adaptable=false;
}

Network::Network(int netid, bool adaptval) {
	name=0; //Defaults to no name                                                                                               
	numnodes=-1;
	numlinks=-1;
	net_id=netid;
	adaptable=adaptval;
}


Network::Network(const Network& network)
{
	std::vector<NNode*>::const_iterator curnode;

	// Copy all the inputs
	for(curnode = network.inputs.begin(); curnode != network.inputs.end(); ++curnode) {
		NNode* n = new NNode(**curnode);
		inputs.push_back(n);
		all_nodes.push_back(n);
	}

	// Copy all the outputs
	for(curnode = network.outputs.begin(); curnode != network.outputs.end(); ++curnode) {
		NNode* n = new NNode(**curnode);
		outputs.push_back(n);
		all_nodes.push_back(n);
	}

	if(network.name)
		name = strdup(network.name);
	else
		name = 0;

	numnodes = network.numnodes;
	numlinks = network.numlinks;
	net_id = network.net_id;
	adaptable = network.adaptable;
}

Network::~Network() {
	if (name!=0)
		delete [] name;

	destroy();  // Kill off all the nodes and links

}

// Puts the network back into an initial state
void Network::flush() {
	std::vector<NNode*>::iterator curnode;

	for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {
		(*curnode)->flushback();
	}
}

// Debugger: Checks network state
void Network::flush_check() {
	std::vector<NNode*>::iterator curnode;
	std::vector<NNode*>::iterator location;
	std::vector<NNode*> seenlist;  //List of nodes not to doublecount

	for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {    
		location= std::find(seenlist.begin(),seenlist.end(),(*curnode));
		if (location==seenlist.end()) {
			seenlist.push_back(*curnode);
			(*curnode)->flushback_check(seenlist);
		}
	}
}

// If all output are not active then return true
bool Network::outputsoff() {
	std::vector<NNode*>::iterator curnode;

	for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {
		if (((*curnode)->activation_count)==0) return true;
	}

	return false;
}

// Print the connections weights to a file separated by only carriage returns
void Network::print_links_tofile(char *filename) {
	std::vector<NNode*>::iterator curnode;
	std::vector<Link*>::iterator curlink;

	std::ofstream oFile(filename);

	//Make sure it worked
	//if (!oFile) {
	//	cerr<<"Can't open "<<filename<<" for output"<<endl;
	//return 0;
	//}

	for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {
		if (((*curnode)->type)!=SENSOR) {
			for(curlink=((*curnode)->incoming).begin(); curlink!=((*curnode)->incoming).end(); ++curlink) {
				oFile << (*curlink)->in_node->node_id << " -> " <<( *curlink)->out_node->node_id << " : " << (*curlink)->weight << std::endl;
			} // end for loop on links
		} //end if
	} //end for loop on nodes

	oFile.close();

} //print_links_tofile

bool Network::init_ctrnn() {
	std::vector<NNode*>::iterator curnode;
	step = 0;
	flush();
	for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {
		if (((*curnode)->type)!=SENSOR) {
			(*curnode)->output=NEAT::fsigmoid((*curnode)->bias,1.0,1.0);
		}

	}

	return true;
}

bool Network::activate_ctrnn(double dt) {
	std::vector<NNode*>::iterator curnode;
	std::vector<Link*>::iterator curlink;
	double add_amount;  //For adding to the activesum

	step++;

	// For each node, compute the sum of its incoming activation
	for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {
		//Ignore SENSORS

		//cout<<"On node "<<(*curnode)->node_id<<endl;

		if (((*curnode)->type)!=SENSOR) {
			(*curnode)->activesum=0;
			(*curnode)->modsum=0;

			// For each incoming connection, add the activity from the connection to the activesum
			for(curlink=((*curnode)->incoming).begin();curlink!=((*curnode)->incoming).end();++curlink) {
				//Handle possible time delays

				if  ( (*curlink)->enabled)
				{
					add_amount=((*curlink)->weight)*(((*curlink)->in_node)->get_output());

					//if ( (*curlink)->mod<= 0.0)
						(*curnode)->activesum+=add_amount;
					//else
					//	(*curnode)->modsum+=add_amount;
				}

			} //End for over incoming links

		} //End if (((*curnode)->type)!=SENSOR)

	} //End for over all nodes
	double delta;
	// Now activate all the non-sensor nodes off their incoming activation
	for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {

		if (((*curnode)->type)!=SENSOR && (*curnode)->incoming.size()!=0) {
			//Only activate if some active input came in
			//cout<<"Activating "<<(*curnode)->node_id<<" with "<<(*curnode)->activesum<<": ";

			//Keep a memory of activations for potential time delayed connections
			(*curnode)->last_activation2=(*curnode)->last_activation;
			(*curnode)->last_activation=(*curnode)->activation;

			//std::cout<<(*curnode)->last_activation<<"\t";
			//Now run the net activation through an activation function
			//if ((*curnode)->ftype==SIGMOID)
			//	(*curnode)->activation=NEAT::fsigmoid((*curnode)->activesum,4.924273,2.4621365);  //Sigmoidal activation- see comments under fsigmoid
			delta = dt/(*curnode)->time_const * (-(*curnode)->activation + (*curnode)->activesum);
			(*curnode)->activation+=delta;

			


			//TOOD between 0..1 or -1 or 1
	//!!!!		(*curnode)->mod_activation=std::tanh( (*curnode)->bias+(*curnode)->modsum);// NEAT::fsigmoid( (*curnode)->bias+(*curnode)->modsum, 1.0, 1.0);;
			//std::cout << randfloat()<<"\t";

			//TODO CHECK IF scaling works correctly
			(*curnode)->output =NEAT::fsigmoid( (*curnode)->bias+(*curnode)->activation, 1.0, 1.0);
			if ((*curnode)->gen_node_label == OUTPUT)
			{
				//(*fnetout) <<(*curnode)->output<<"\t";
			}

			//if ( (*curnode)->bias <= 0.0f && (*curnode)->gen_node_label == HIDDEN)
			//{
			//	//REMOVE!!!!!!!!!!!!!!
			//	(*curnode)->output = 0.0f;
			//}

			//Step might not be the same at the start???

		//	(*curnode)->output =NEAT::fsigmoid( sin(step/10.0+ (*curnode)->bias )*3.2+(*curnode)->activation, 1.0, 1.0);
			
			
			//TODO ******************tanh((*curnode)->bias+(*curnode)->activation);

			//if (NEAT::HYPERNEAT_SETTING==4) //SATURATE
			
			
	//	(*curnode)->output += randfloat()*0.2-0.1;

			if ((*curnode)->output<0.0)
			{
			//	std::cout<<"output negative";
				(*curnode)->output = 0.0;
			}

				//NEAT::fsigmoid( (*curnode)->bias+(*curnode)->activation, 1.0, 1.0)+ (randfloat()*0.2-0.1); 
			// tanh( (*curnode)->bias+(*curnode)->activation);//+ (randfloat()*0.2-0.1);
			//NEAT::fsigmoid((*curnode)->bias+(*curnode)->activation + (randfloat()*0.2-0.1),1.0,1.0);// + (randfloat()*0.2-0.1);
			//tanh((*curnode)->bias+(*curnode)->activation) + (randfloat()*0.2-0.1);/
			//RISI
			//cout<<(*curnode)->activation<<endl;

			//Increment the activation_count
			//First activation cannot be from nothing!!
			(*curnode)->activation_count++;

		}
	}
	float oldWeight, weight_sign;
	//(*fnetout)<<"\n";

	float mod_signal = 1.0f; //TODO give negative reward with mouseclick

	return true;


	/*static int timer=0;
	timer++;
	if (timer>350)
	{
		mod_signal= -1.0f;
		std::cout<<"*";
	}*/

	//if (NEAT::HYPERNEAT_SETTING==2 || NEAT::HYPERNEAT_SETTING==4) { //adaptation enabled

		//std::cout << "ADAPTING" << std:endl;
		
	//float add_amount=0.0f;

		/*if ((*( (*all_nodes.begin())->incoming).begin())->trace==0)
			(*( (*all_nodes.begin())->incoming).begin())->trace=1;

		(*( (*all_nodes.begin())->incoming).begin())->trace = 0.995 * (*( (*all_nodes.begin())->incoming).begin())->trace;

		std::cout<<"Trace "<<(*( (*all_nodes.begin())->incoming).begin())->trace<<"\n";*/

		//return true;

	float B = 0.995f; 

		// ADAPTATION:  Adapt weights based on activations 
		for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {
			//Ignore SENSORS

			//cout<<"On node "<<(*curnode)->node_id<<endl;
			// For each incoming connection, perform adaptation based on the trait of the connection 
			for(curlink=((*curnode)->incoming).begin();curlink!=((*curnode)->incoming).end();++curlink) {

				//if ( !(*curlink)->crossmodule) continue;

				//if (((*curnode)->type)!=SENSOR) {
	//!			if ( (*curnode)->gen_node_label==OUTPUT) {  //Only adapt weights going to the outputs so we don't mess with the occilatory pattern of the hidden nodes
				//	std::cout<<(*curlink)->A<<"\t"<<(*cur link)->B<<"\t"<<(*curlink)->C<<"\t"<<(*curlink)->lrate<<"\n";
					
				//	std::cout<<(*curlink)->weight<<"\t";

					//old adaptation
					if (NEAT::HYPERNEAT_SETTING==2) //ABC
					{
						if ( (*curlink)->mod <= 0) //normal connection
						{
							
							(*curlink)->weight+=(*curlink)->out_node->mod_activation * 
								(*curlink)->lrate*5.0f*
								((*curlink)->A* (*curlink)->in_node->get_active_out()* (*curlink)->out_node->get_active_out()
								+ (*curlink)->B*(*curlink)->in_node->get_active_out()
								+ (*curlink)->C* (*curlink)->out_node->get_active_out() + (*curlink)->D);
						}
					}
					//Cacluate trace
					//get_active_out()
				//	std::cout<<(*curlink)->in_node->get_active_out()<<"\t"<<(*curlink)->out_node->output<<"\n";
					add_amount = 0.0f;
					//std::cout<<( (*curlink)->in_node->output * (*curlink)->out_node->output )<<"\n";
					if ((*curlink)->in_node->output * (*curlink)->out_node->output > 1.2)
					{
						add_amount = 0.5f;
					}
					else if ( (*curlink)->in_node->output * (*curlink)->out_node->output < 0.00001) 
					{
							add_amount = -1.0;
							//std::cout<<( (*curlink)->in_node->output * (*curlink)->out_node->output )<<"\n";
					
					}

					(*curlink)->trace = B * (*curlink)->trace + add_amount;//(1.0f - B) * add_amount;
					if ((*curlink)->trace< 0.0) (*curlink)->trace = 0.0;

					//TODO maybe make activation function non bimodal
					if (NEAT::HYPERNEAT_SETTING==4) //SATURATE
					{
						weight_sign = 1.0f;
						if ((*curlink)->weight<0.0) weight_sign=-1.0f;

						oldWeight = (*curlink)->weight;
						(*curlink)->weight+=weight_sign * mod_signal * (*curlink)->in_node->get_active_out() * (*curlink)->out_node->get_active_out() + (randfloat()*0.2-0.1);
					
						//if (curnode==all_nodes.begin())
						//	std::cout<<oldWeight<<"\t"<<weight_sign<<"\n";
				
						if (oldWeight>0 && (*curlink)->weight < 0) (*curlink)->weight = (randfloat()*0.1);//0.01;
						if (oldWeight<0 && (*curlink)->weight > 0) (*curlink)->weight = -(randfloat()*0.1);//-0.01;
					}
					
					if ((*curlink)->weight > 10) (*curlink)->weight=10;
					if ((*curlink)->weight < -10) (*curlink)->weight=-10;

					//std::cout<<(*curlink)->weight<<"\n";

					//		      hebbian((*curlink)->weight,maxweight,
					//			      (*curlink)->in_node->get_active_out(), 
					//			      (*curlink)->out_node->get_active_out(),
	//!			}

			}
	//		std::cout<<"\n";
		}

	//}
	return true;
}

void Network::resetCTRNN()
{
	//init_ctrnn();
	//for(int x=0;x<50;x++)
	//		activate_ctrnn(0.02);

	std::vector<NNode*>::iterator curnode;
	std::vector<Link*>::iterator curlink;
	for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {

		//	(*curnode)->time_backup = (*curnode)->time_const;
	//		(*curnode)->bias_backup = (*curnode)->bias;
				//NEW
			(*curnode)->activation = 0.0;
			(*curnode)->last_activation = 0.0;
	}
}

void Network::revertChanges()
{
		std::vector<NNode*>::iterator curnode;
	std::vector<Link*>::iterator curlink;
	for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {

			(*curnode)->time_const = (*curnode)->time_backup;
			(*curnode)->bias = (*curnode)->bias_backup;

			
	}
}

void Network::reinforce(double r)
{
	//return;
	std::vector<NNode*>::iterator curnode;
	std::vector<Link*>::iterator curlink;
	for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {

			(*curnode)->time_backup = (*curnode)->time_const;
			(*curnode)->bias_backup = (*curnode)->bias;
			//	//NEW
			//(*curnode)->activation = 0.0;
			//(*curnode)->last_activation = 0.0;

			if ( (*curnode)->gen_node_label==HIDDEN)
			{
				(*curnode)->time_const += randfloat()*r-r/2.0;
				(*curnode)->bias +=  randfloat()*r-r/2.0;
				if ((*curnode)->time_const<0.1) (*curnode)->time_const = 0.1;
				if ((*curnode)->time_const>2.0) (*curnode)->time_const= 2.0;

				if ((*curnode)->bias<-3.0) (*curnode)->bias = -3.0;
				if ((*curnode)->bias>3.0) (*curnode)->bias= 3.0;
				//ALSO BIAS?
			}

		//for(curlink=((*curnode)->incoming).begin();curlink!=((*curnode)->incoming).end();++curlink) {

		//	if ((*curlink)->lrate!=0.0)
		//	{
		//		(*curlink)->weight+= (*curlink)->trace * r * 
		//							(*curlink)->lrate*1.0f* (*curlink)->in_node->get_active_out()* (*curlink)->out_node->get_active_out();

		//				//(*curlink)->weight+= (*curlink)->trace * r * 
		//				//			(*curlink)->lrate*3.0f*
		//				//			((*curlink)->A* (*curlink)->in_node->get_active_out()* (*curlink)->out_node->get_active_out()
		//				//			+ (*curlink)->B*(*curlink)->in_node->get_active_out()
		//				//			+ (*curlink)->C* (*curlink)->out_node->get_active_out() + (*curlink)->D);

		//		if ( (*curlink)->weight > 5.0) (*curlink)->weight = 5.0;
		//		if ( (*curlink)->weight < -5.0) (*curlink)->weight = -5.0;
		//	}
		//	
		//	//if ( !(*curlink)->crossmodule ) continue;

		//	//(*curlink)->weight -= 0.3 * (*curlink)->weight * (*curlink)->trace;

		//	//if ( (*curlink)->weight >= 0.0)
		//	//	(*curlink)->weight +=r * (*curlink)->trace;
		//	//else
		//	//	(*curlink)->weight -=r * (*curlink)->trace;

		//	//Weights only increase?

		//	//std::cout<<"R = "<< ( r * (*curlink)->trace)<<"\n";
		//	//(*curlink)->weight += 0.3 * (*curlink)->trace*r; //A was 0.3
		//	//if ( (*curlink)->trace < 0.0)
		//	//	std::cout <<" "<<(*curlink)->trace<<"\n";
		//}
	}
}

// Activates the net such that all outputs are active
// Returns true on success;
bool Network::activate() {
	std::vector<NNode*>::iterator curnode;
	std::vector<Link*>::iterator curlink;
	//double add_amount;  //For adding to the activesum
	//bool onetime; //Make sure we at least activate once
	abortcount=0;  //Used in case the output is somehow truncated from the network

	//cout<<"Activating network: "<<this->genotype<<endl;

	//Keep activating until all the outputs have become active 
	//(This only happens on the first activation, because after that they
	// are always active)

	onetime=false;

	while(outputsoff()||!onetime) {

		++abortcount;

		if (abortcount==20) {
			return false;
			//cout<<"Inputs disconnected from output!"<<endl;
		}
		//std::cout<<"Outputs are off"<<std::endl;

		// For each node, compute the sum of its incoming activation 
		for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {
			//Ignore SENSORS

			//cout<<"On node "<<(*curnode)->node_id<<endl;

			if (((*curnode)->type)!=SENSOR) {
				(*curnode)->activesum=0;
				(*curnode)->active_flag=false;  //This will tell us if it has any active inputs

				// For each incoming connection, add the activity from the connection to the activesum 
				for(curlink=((*curnode)->incoming).begin();curlink!=((*curnode)->incoming).end();++curlink) {
					//Handle possible time delays
					if (!((*curlink)->time_delay)) {
						add_amount=((*curlink)->weight)*(((*curlink)->in_node)->get_active_out());
						if ((((*curlink)->in_node)->active_flag)||
							(((*curlink)->in_node)->type==SENSOR)) (*curnode)->active_flag=true;
						(*curnode)->activesum+=add_amount;
						//std::cout<<"Node "<<(*curnode)->node_id<<" adding "<<add_amount<<" from node "<<((*curlink)->in_node)->node_id<<std::endl;
					}
					else {
						//Input over a time delayed connection
						add_amount=((*curlink)->weight)*(((*curlink)->in_node)->get_active_out_td());
						(*curnode)->activesum+=add_amount;
					}

				} //End for over incoming links

			} //End if (((*curnode)->type)!=SENSOR) 

		} //End for over all nodes

		// Now activate all the non-sensor nodes off their incoming activation 
		for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {

			if (((*curnode)->type)!=SENSOR) {
				//Only activate if some active input came in
				if ((*curnode)->active_flag) {
					//cout<<"Activating "<<(*curnode)->node_id<<" with "<<(*curnode)->activesum<<": ";

					//Keep a memory of activations for potential time delayed connections
					(*curnode)->last_activation2=(*curnode)->last_activation;
					(*curnode)->last_activation=(*curnode)->activation;

					//If the node is being overriden from outside,
					//stick in the override value
					if ((*curnode)->overridden()) {
						//Set activation to the override value and turn off override
						(*curnode)->activate_override();
					}
					else {
						//Now run the net activation through an activation function
						switch((*curnode)->ftype)
						{
						case ACTIVATION_FUNCTION_TANH:
							//Sigmoidal activation- see comments under fsigmoid
							(*curnode)->activation=tanh((*curnode)->activesum); //NEAT::fsigmoid((*curnode)->activesum,4.924273,2.4621365);
							break;
						case ACTIVATION_FUNCTION_SIGMOID:
							//Sigmoidal activation- see comments under fsigmoid
							(*curnode)->activation=fsigmoid((*curnode)->activesum,4.924273,2.4621365);
							break;
						case ACTIVATION_FUNCTION_SIN:
							(*curnode)->activation=NEAT::fsin((*curnode)->activesum);
							break;
						case ACTIVATION_FUNCTION_COS:
							(*curnode)->activation=NEAT::fcos((*curnode)->activesum);
							break;
						case ACTIVATION_FUNCTION_GAUSSIAN:
							(*curnode)->activation=NEAT::fgaussian((*curnode)->activesum);
							break;
						case ACTIVATION_FUNCTION_ABS:
							(*curnode)->activation=abs((*curnode)->activesum);
							break;
						case ACTIVATION_FUNCTION_RECT:
							if ( abs ( (*curnode)->activesum) <0.5) (*curnode)->activation=1;
							else (*curnode)->activation=0;
							break;
						//case ACTIVATION_FUNCTION_SQUARE:
						//	(*curnode)->activation=NEAT::fsquare((*curnode)->activesum);
						//	break;
						//case ACTIVATION_FUNCTION_ABS_ROOT:
						//	(*curnode*)->activation=NEAT::fabsroot((*curnode)->activesum);
						//	break;
					//	case ACTIVATION_FUNCTION_LINEAR:
						//	(*curnode)->activation=NEAT::flinear((*curnode)->activesum);
						//	break;
						}
					}
					//cout<<(*curnode)->activation<<endl;

					//Increment the activation_count
					//First activation cannot be from nothing!!
					(*curnode)->activation_count++;
				}
			}
		}

		onetime=true;
	}

	if (false && adaptable) {

		//std::cout << "ADAPTING" << std:endl;

		// ADAPTATION:  Adapt weights based on activations 
		for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {
			//Ignore SENSORS

			//cout<<"On node "<<(*curnode)->node_id<<endl;

			if (((*curnode)->type)!=SENSOR) {

				// For each incoming connection, perform adaptation based on the trait of the connection 
				for(curlink=((*curnode)->incoming).begin();curlink!=((*curnode)->incoming).end();++curlink) {

					if (((*curlink)->trait_id==2)||
						((*curlink)->trait_id==3)||
						((*curlink)->trait_id==4)) {

							//In the recurrent case we must take the last activation of the input for calculating hebbian changes
							if ((*curlink)->is_recurrent) {
								(*curlink)->weight=
									hebbian((*curlink)->weight,maxweight,
									(*curlink)->in_node->last_activation, 
									(*curlink)->out_node->get_active_out(),
									(*curlink)->params[0],(*curlink)->params[1],
									(*curlink)->params[2]);


							}
							else { //non-recurrent case
								(*curlink)->weight=
									hebbian((*curlink)->weight,maxweight,
									(*curlink)->in_node->get_active_out(), 
									(*curlink)->out_node->get_active_out(),
									(*curlink)->params[0],(*curlink)->params[1],
									(*curlink)->params[2]);
							}
					}

				}

			}

		}

	} //end if (adaptable)

	return true;  
}
void Network::print_to_SVG(const char* filename)
{
	int xdisplace = 1400;
	int ydisplace = 800;
	float dtx = 110; //was 1100
	float dty = 100; //was 1000
	int bigcircle = 16;
	int smallcircle = 10;

	std::ofstream fout(filename);

	//saveAsEPS = true;
	float width;
	//    String fill1, fill2;
	//    String epsStroke;
	char* epsStroke = "\"black\"";

	fout<<"<?xml version=\"1.0\" standalone=\"no\"?>\n";
	fout<<"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
	fout<<"<svg width=\"12cm\" height=\"12cm\" viewBox=\"0 0 3200 2600\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n";
	fout<<"<desc>Network differences</desc>\n";
	fout<<"<defs>\n";
	fout<<"<marker id=\"Triangle\" viewBox=\"0 0 10 10\" refX=\"0\" refY=\"5\" markerUnits=\"strokeWidth\" markerWidth=\"4\" markerHeight=\"3\" orient=\"auto\"> <path d=\"M 0 0 L 10 5 L 0 10 z\" />\n";
	fout<<"</marker>\n";
	fout<<"</defs>\n";

	std::vector<NNode*>::iterator curnode;
	std::vector<Link*>::iterator curlink;

	float x, y;
	float x1, y1, x2, y2;
	//Inter module connections first
	for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) 
	{
		x = (*curnode)->x_sub_pos+(*curnode)->x_global*3.0;
		y = (*curnode)->y_sub_pos+(*curnode)->y_global*3.0;

		// For each incoming connection, add the activity from the connection to the activesum 
		for(curlink=((*curnode)->incoming).begin();curlink!=((*curnode)->incoming).end();++curlink) 
		{
			if ( !(*curlink)->enabled) 
			{
				std::cout<<"Connection disabled "<<(*curlink)->in_node->node_id<<" "<<(*curlink)->out_node->node_id<<"\n";
		
				continue;
			}

			float width= abs((*curlink)->weight);

			if ( (*curlink)->in_node->x_global != (*curlink)->out_node->x_global
				||  (*curlink)->in_node->y_global != (*curlink)->out_node->y_global)
			{
				epsStroke = "\"green\"";
			//	width = 10;
				
			} else continue;
			
			x1 = (*curlink)->in_node->x_sub_pos+(*curlink)->in_node->x_global*3.0;
			y1 = (*curlink)->in_node->y_sub_pos+ (*curlink)->in_node->y_global*3.0;

			x2 = (*curlink)->out_node->x_sub_pos+(*curlink)->out_node->x_global*3.0;
			y2 = (*curlink)->out_node->y_sub_pos+ (*curlink)->out_node->y_global*3.0;

			fout<<"<line x1=\"" << ( x1 * dtx + xdisplace) << "\" y1=\"" << (ydisplace - y1 * dty) <<
				"\" x2=\"" << (x2 * dtx + xdisplace) << "\" y2=\"" << (ydisplace - y2 * dty) << "\" stroke-width=\"" << width <<
				"\"  stroke=" << epsStroke << "/>\n";
		}
				
		//for(curlink=((*curnode)->outgoing).begin();curlink!=((*curnode)->outgoing).end();++curlink) 
		//{
		//	if ( !(*curlink)->enabled) continue;

		//	float width= (*curlink)->weight;

		//	if ( (*curlink)->in_node->x_global != (*curlink)->out_node->x_global
		//		||  (*curlink)->in_node->y_global != (*curlink)->out_node->y_global)
		//	{
		//		epsStroke = "\"green\"";
		//		
		//		width = 10;
		//		
		//	} else continue;

		//	x1 = (*curlink)->in_node->x_sub_pos+(*curlink)->in_node->x_global*2.0;
		//	y1 = (*curlink)->in_node->y_sub_pos+ (*curlink)->in_node->y_global*2.0;

		//	x2 = (*curlink)->out_node->x_sub_pos+(*curlink)->out_node->x_global*2.0;
		//	y2 = (*curlink)->out_node->y_sub_pos+ (*curlink)->out_node->y_global*2.0;

		//	//fout<<"<line x1=\"" << ( (*curlink)->in_node->x_sub_pos * dtx + xdisplace) << "\" y1=\"" << (ydisplace - (*curlink)->in_node->y_sub_pos * dty) <<
		//	//	"\" x2=\"" << ( (*curlink)->out_node->x_sub_pos * dtx + xdisplace) << "\" y2=\"" << (ydisplace - (*curlink)->out_node->y_sub_pos * dty) << "\" stroke-width=\"" << width <<
		//	//	"\"  stroke=" << epsStroke << "/>\n";

		//	fout<<"<line x1=\"" << ( x1 * dtx + xdisplace) << "\" y1=\"" << (ydisplace - y1 * dty) <<
		//		"\" x2=\"" << (x2 * dtx + xdisplace) << "\" y2=\"" << (ydisplace - y2 * dty) << "\" stroke-width=\"" << width <<
		//		"\"  stroke=" << epsStroke << "/>\n";

		//}
	}

		//Now all the other connections
	for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) 
	{
	//	if ( (*curnode)->bias<=0.0 && (*curnode)->gen_node_label == HIDDEN) continue;

		//std::cout<<(*curnode)->x_global<<"\n"<<(*curnode)->y_global<<"\n";

		x = (*curnode)->x_sub_pos+(*curnode)->x_global*3.0;
		y = (*curnode)->y_sub_pos+(*curnode)->y_global*3.0;

		// For each incoming connection, add the activity from the connection to the activesum 
		for(curlink=((*curnode)->incoming).begin();curlink!=((*curnode)->incoming).end();++curlink) 
		{
			if ( !(*curlink)->enabled) 
			{
				std::cout<<"Connection disabled "<<(*curlink)->in_node->node_id<<" "<<(*curlink)->out_node->node_id<<"\n";
		
				continue;
			}

			float width= abs((*curlink)->weight);

			if (width >= 0.0f)
			{
				//penConnection.Color = Color.Black;
				epsStroke = "\"black\"";
			}
			else
			{
				// penConnection.Color = Color.Red;
				epsStroke = "\"red\"";
			}

			if ( (*curlink)->in_node->x_global != (*curlink)->out_node->x_global
				||  (*curlink)->in_node->y_global != (*curlink)->out_node->y_global)
			{
				epsStroke = "\"green\"";
				continue;
				width = 10;
				
			}

			x1 = (*curlink)->in_node->x_sub_pos+(*curlink)->in_node->x_global*3.0;
			y1 = (*curlink)->in_node->y_sub_pos+ (*curlink)->in_node->y_global*3.0;

			x2 = (*curlink)->out_node->x_sub_pos+(*curlink)->out_node->x_global*3.0;
			y2 = (*curlink)->out_node->y_sub_pos+ (*curlink)->out_node->y_global*3.0;

			fout<<"<line x1=\"" << ( x1 * dtx + xdisplace) << "\" y1=\"" << (ydisplace - y1 * dty) <<
				"\" x2=\"" << (x2 * dtx + xdisplace) << "\" y2=\"" << (ydisplace - y2 * dty) << "\" stroke-width=\"" << width <<
				"\"  stroke=" << epsStroke << "/>\n";

		}

		for(curlink=((*curnode)->outgoing).begin();curlink!=((*curnode)->outgoing).end();++curlink) 
		{
			if ( !(*curlink)->enabled) continue;

			float width= (*curlink)->weight;

			if (width >= 0.0f)
			{
				//penConnection.Color = Color.Black;
				epsStroke = "\"black\"";
			}
			else
			{
				// penConnection.Color = Color.Red;
				epsStroke = "\"red\"";
			}
			if ( (*curlink)->in_node->x_global != (*curlink)->out_node->x_global
				||  (*curlink)->in_node->y_global != (*curlink)->out_node->y_global)
			{
				epsStroke = "\"green\"";
				continue;
				width = 10;
				
			}
			x1 = (*curlink)->in_node->x_sub_pos+(*curlink)->in_node->x_global*3.0;
			y1 = (*curlink)->in_node->y_sub_pos+ (*curlink)->in_node->y_global*3.0;

			x2 = (*curlink)->out_node->x_sub_pos+(*curlink)->out_node->x_global*3.0;
			y2 = (*curlink)->out_node->y_sub_pos+ (*curlink)->out_node->y_global*3.0;

			//fout<<"<line x1=\"" << ( (*curlink)->in_node->x_sub_pos * dtx + xdisplace) << "\" y1=\"" << (ydisplace - (*curlink)->in_node->y_sub_pos * dty) <<
			//	"\" x2=\"" << ( (*curlink)->out_node->x_sub_pos * dtx + xdisplace) << "\" y2=\"" << (ydisplace - (*curlink)->out_node->y_sub_pos * dty) << "\" stroke-width=\"" << width <<
			//	"\"  stroke=" << epsStroke << "/>\n";

			fout<<"<line x1=\"" << ( x1 * dtx + xdisplace) << "\" y1=\"" << (ydisplace - y1 * dty) <<
				"\" x2=\"" << (x2 * dtx + xdisplace) << "\" y2=\"" << (ydisplace - y2 * dty) << "\" stroke-width=\"" << width <<
				"\"  stroke=" << epsStroke << "/>\n";

		}
	}
		
	//Now the nodes on top
	for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) 
		{
		//	if ( (*curnode)->bias<=0.0 && (*curnode)->gen_node_label == HIDDEN) continue;

			//std::cout<<(*curnode)->x_global<<"\n"<<(*curnode)->y_global<<"\n";

			x = (*curnode)->x_sub_pos+(*curnode)->x_global*3.0;
			y = (*curnode)->y_sub_pos+(*curnode)->y_global*3.0;


			if ( (*curnode)->gen_node_label==INPUT)
			{
				fout<<"<circle cx=\"" << ( x* dtx + xdisplace) << "\" cy=\"" << (ydisplace - y * dty) << "\" r=\"25\" stroke-width=\"3\" stroke=" << "\"black\"" << " fill=\"" << "black" << "\"/>\n";
			}
			else if ( (*curnode)->gen_node_label==OUTPUT)
			{
				fout<<"<circle cx=\"" << ( x * dtx + xdisplace) << "\" cy=\"" << (ydisplace - y * dty) << "\" r=\"25\" stroke-width=\"3\" stroke=" << "\"black\"" << " fill=\"" << "red" << "\"/>\n";
			}
			else
			{
				fout<<"<circle cx=\"" << ( x * dtx + xdisplace) << "\" cy=\"" << (ydisplace - y * dty) << "\" r=\"25\" stroke-width=\"3\" stroke=" << "\"black\"" << " fill=\"" << "yellow" << "\"/>\n";
			}
		}

	//SW.WriteLine("<rect x=\"1\" y=\"1\" width=\"1198\" height=\"398\" fill=\"none\" stroke=\"blue\" stroke-width=\"2\"/>");
	//SW.WriteLine("<rect x=\"400\" y=\"100\" width=\"400\" height=\"200\" fill=\"yellow\" stroke=\"navy\" stroke-width=\"10\"/>");
	fout<<"</svg>"<<std::endl;
	// saveAsEPS = false;
	fout.close();
}

// THIS WAS NOT USED IN THE FINAL VERSION, AND NOT FULLY IMPLEMENTED,   
// BUT IT SHOWS HOW SOMETHING LIKE THIS COULD BE INITIATED
// Note that checking networks for loops in general in not necessary
// and therefore I stopped writing this function
// Check Network for loops.  Return true if its ok, false if there is a loop.
//bool Network::integrity() {
//  std::vector<NNode*>::iterator curnode;
//  std::vector<std::vector<NNode*>*> paths;
//  int count;
//  std::vector<NNode*> *newpath;
//  std::vector<std::vector<NNode*>*>::iterator curpath;

//  for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {
//    newpath=new std::vector<NNode*>();
//    paths.push_back(newpath);
//    if (!((*curnode)->integrity(newpath))) return false;
//  }

//Delete the paths now that we are done
//  curpath=paths.begin();
//  for(count=0;count<paths.size();count++) {
//    delete (*curpath);
//    curpath++;
//  }

//  return true;
//}

// Prints the values of its outputs
void Network::show_output() {
	std::vector<NNode*>::iterator curnode;
	int count;

	//if (name!=0)
	//  cout<<"Network "<<name<<" with id "<<net_id<<" outputs: (";
	//else cout<<"Network id "<<net_id<<" outputs: (";

	count=1;
	for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {
		std::cout<<"[Output #"<<count<<": "<<(*curnode)->output<< " " << (*curnode)->time_const << " " << (*curnode)->bias << "] " << std::endl;
		count++;
	}

}

void Network::show_activation() {
	std::vector<NNode*>::iterator curnode;
	int count;

	//if (name!=0)
	//  cout<<"Network "<<name<<" with id "<<net_id<<" outputs: (";
	//else cout<<"Network id "<<net_id<<" outputs: (";

	count=1;
	for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {
		//cout<<"[Output #"<<count<<": "<<(*curnode)<<"] ";
		count++;
	}

	//cout<<")"<<endl;
}

void Network::show_input() {
	std::vector<NNode*>::iterator curnode;
	int count;

	//if (name!=0)
	//  cout<<"Network "<<name<<" with id "<<net_id<<" inputs: (";
	//else cout<<"Network id "<<net_id<<" outputs: (";

	count=1;
	for(curnode=inputs.begin();curnode!=inputs.end();++curnode) {
		//cout<<"[Input #"<<count<<": "<<(*curnode)<<"] ";
		count++;
	}

	//cout<<")"<<endl;
}

// Add an input
void Network::add_input(NNode *in_node) {
	inputs.push_back(in_node);
}

// Add an output
void Network::add_output(NNode *out_node) {
	outputs.push_back(out_node);
}

// Takes an array of sensor values and loads it into SENSOR inputs ONLY
void Network::load_sensors(double *sensvals) {
	//int counter=0;  //counter to move through array
	std::vector<NNode*>::iterator sensPtr;

	for(sensPtr=inputs.begin();sensPtr!=inputs.end();++sensPtr) {
		//only load values into SENSORS (not BIASes)
		if (((*sensPtr)->type)==SENSOR) {
			(*sensPtr)->sensor_load(*sensvals);
			sensvals++;
		}
	}
}

void Network::load_sensors(const std::vector<float> &sensvals) {
	//int counter=0;  //counter to move through array
	std::vector<NNode*>::iterator sensPtr;
	std::vector<float>::const_iterator valPtr;

	for(valPtr = sensvals.begin(), sensPtr = inputs.begin(); sensPtr != inputs.end() && valPtr != sensvals.end(); ++sensPtr, ++valPtr) {
		//only load values into SENSORS (not BIASes)
		if (((*sensPtr)->type)==SENSOR) {
			(*sensPtr)->sensor_load(*valPtr);
			//sensvals++;
		}
	}
}


// Takes and array of output activations and OVERRIDES 
// the outputs' actual activations with these values (for adaptation)
void Network::override_outputs(double* outvals) {

	std::vector<NNode*>::iterator outPtr;

	for(outPtr=outputs.begin();outPtr!=outputs.end();++outPtr) {
		(*outPtr)->override_output(*outvals);
		outvals++;
	}

}

void Network::give_name(char *newname) {
	char *temp;
	char *temp2;
	temp=new char[strlen(newname)+1];
	strcpy(temp,newname);
	if (name==0) name=temp;
	else {
		temp2=name;
		delete temp2;
		name=temp;
	}
}

// The following two methods recurse through a network from outputs
// down in order to count the number of nodes and links in the network.
// This can be useful for debugging genotype->phenotype spawning 
// (to make sure their counts correspond)

int Network::nodecount() {
	int counter=0;
	std::vector<NNode*>::iterator curnode;
	std::vector<NNode*>::iterator location;
	std::vector<NNode*> seenlist;  //List of nodes not to doublecount

	for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {

		location = std::find(seenlist.begin(),seenlist.end(),(*curnode));
		if (location==seenlist.end()) {
			counter++;
			seenlist.push_back(*curnode);
			nodecounthelper((*curnode),counter,seenlist);
		}
	}

	numnodes=counter;

	return counter;

}

void Network::nodecounthelper(NNode *curnode,int &counter,std::vector<NNode*> &seenlist) {
	std::vector<Link*> innodes=curnode->incoming;
	std::vector<Link*>::iterator curlink;
	std::vector<NNode*>::iterator location;

	if (!((curnode->type)==SENSOR)) {
		for(curlink=innodes.begin();curlink!=innodes.end();++curlink) {
			location= std::find(seenlist.begin(),seenlist.end(),((*curlink)->in_node));
			if (location==seenlist.end()) {
				counter++;
				seenlist.push_back((*curlink)->in_node);
				nodecounthelper((*curlink)->in_node,counter,seenlist);
			}
		}

	}

}

int Network::linkcount() {
	int counter=0;
	std::vector<NNode*>::iterator curnode;
	std::vector<NNode*> seenlist;  //List of nodes not to doublecount

	for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {
		linkcounthelper((*curnode),counter,seenlist);
	}

	numlinks=counter;

	return counter;

}

void Network::linkcounthelper(NNode *curnode,int &counter,std::vector<NNode*> &seenlist) {
	std::vector<Link*> inlinks=curnode->incoming;
	std::vector<Link*>::iterator curlink;
	std::vector<NNode*>::iterator location;

	location = std::find(seenlist.begin(),seenlist.end(),curnode);
	if ((!((curnode->type)==SENSOR))&&(location==seenlist.end())) {
		seenlist.push_back(curnode);

		for(curlink=inlinks.begin();curlink!=inlinks.end();++curlink) {
			counter++;
			linkcounthelper((*curlink)->in_node,counter,seenlist);
		}

	}

}

// Destroy will find every node in the network and subsequently
// delete them one by one.  Since deleting a node deletes its incoming
// links, all nodes and links associated with a network will be destructed
// Note: Traits are parts of genomes and not networks, so they are not
//       deleted here
void Network::destroy() {
	std::vector<NNode*>::iterator curnode;
	std::vector<NNode*>::iterator location;
	std::vector<NNode*> seenlist;  //List of nodes not to doublecount

	// Erase all nodes from all_nodes list 

	for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {
		delete (*curnode);
	}


	// ----------------------------------- 

	//  OLD WAY-the old way collected the nodes together and then deleted them

	//for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {
	//cout<<seenstd::vector<<endl;
	//cout<<curnode<<endl;
	//cout<<curnode->node_id<<endl;

	//  location=find(seenlist.begin(),seenlist.end(),(*curnode));
	//  if (location==seenlist.end()) {
	//    seenlist.push_back(*curnode);
	//    destroy_helper((*curnode),seenlist);
	//  }
	//}

	//Now destroy the seenlist, which is all the NNodes in the network
	//for(curnode=seenlist.begin();curnode!=seenlist.end();++curnode) {
	//  delete (*curnode);
	//}
}

void Network::destroy_helper(NNode *curnode,std::vector<NNode*> &seenlist) {
	std::vector<Link*> innodes=curnode->incoming;
	std::vector<Link*>::iterator curlink;
	std::vector<NNode*>::iterator location;

	if (!((curnode->type)==SENSOR)) {
		for(curlink=innodes.begin();curlink!=innodes.end();++curlink) {
			location = std::find(seenlist.begin(),seenlist.end(),((*curlink)->in_node));
			if (location==seenlist.end()) {
				seenlist.push_back((*curlink)->in_node);
				destroy_helper((*curlink)->in_node,seenlist);
			}
		}

	}

}

// This checks a POTENTIAL link between a potential in_node and potential out_node to see if it must be recurrent 
bool Network::is_recur(NNode *potin_node,NNode *potout_node,int &count,int thresh) {
	std::vector<Link*>::iterator curlink;


	++count;  //Count the node as visited

	if (count>thresh) {
		//cout<<"returning false"<<endl;
		return false;  //Short out the whole thing- loop detected
	}

	if (potin_node==potout_node) return true;
	else {
		//Check back on all links...
		for(curlink=(potin_node->incoming).begin();curlink!=(potin_node->incoming).end();curlink++) {
			//But skip links that are already recurrent
			//(We want to check back through the forward flow of signals only
			if (!((*curlink)->is_recurrent)) {
				if (is_recur((*curlink)->in_node,potout_node,count,thresh)) return true;
			}
		}
		return false;
	}
}

int Network::input_start() {
	input_iter=inputs.begin();
	return 1;
}

int Network::load_in(double d) {
	(*input_iter)->sensor_load(d);
	input_iter++;
	if (input_iter==inputs.end()) return 0;
	else return 1;
}


//Find the maximum number of neurons between an ouput and an input
int Network::max_depth() {
	std::vector<NNode*>::iterator curoutput; //The current output we are looking at
	int cur_depth; //The depth of the current node
	int max=0; //The max depth

	for(curoutput=outputs.begin();curoutput!=outputs.end();curoutput++) {
		cur_depth=(*curoutput)->depth(0,this);
		if (cur_depth>max) max=cur_depth;
	}

	return max;

}

