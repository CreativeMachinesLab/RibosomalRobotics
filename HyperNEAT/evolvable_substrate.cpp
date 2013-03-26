#include "evolvable_substrate.h"
//#include <hash_set>
#include <iostream>
#include <cmath>

using namespace NEAT;
using namespace std;

//TODO: recursive activation (avoid recursive connections), seeding, leo, adaptation, c
//why not start biped coordinates at -1??? add iterativly from discoverd hidden nodes
//use full hypercube. change in genesis function
////TODO bias for the OUTPUT NODES
//tanh function in node.cpp
//check if scaling works right
//disable repeating loadeing of substrate file


EvolvableSubstrate::EvolvableSubstrate()
{}

ESQuadTree::~ESQuadTree()
{
	std::vector<ESQuadTree*>::iterator curChild;
	for (curChild=childs.begin(); curChild!=childs.end();curChild++)
		delete (*curChild);
}

void ESQuadTree::getConnections(std::vector<Connection2D*> &list, const float &treshold, const float &bandLevel, const float &divisionThreshold, const int &maximumResolution, const int &initialResolution)
{
	float parentVariance = variance(); 

	//cout<<level<<" "<<treshold<<"\n";

	//are we at the level above the leaf nodes and we haven't reached the maximum resolution and the threshold still 
	//higher than threshold then divide further  (level == initialResolution-1) && 
	if (parentVariance > divisionThreshold && (level+1 < maximumResolution) && (childs[0]!=NULL) && (childs[0]->childs.size()==0) )//)// && childs[0]->childs[0]!=NULL )
	{
		//cout << level<<" "<<maximumResolution<<"\n";
		std::vector<ESQuadTree*>::iterator curchild;
		for (curchild=childs.begin();curchild!=childs.end();curchild++)
		{
			(*curchild)->divide( (*curchild)->level + 1);//divide once
		}
	}

	float childVariance;
	float targetx, targety;
	float weight;

	std::vector<ESQuadTree*>::iterator curchild;
	for (curchild=childs.begin();curchild!=childs.end();curchild++)
	{
		childVariance = (*curchild)->variance();

		if (childVariance >= treshold)  //greater or greater and equal?
		{
			(*curchild)->getConnections(list, treshold, bandLevel, divisionThreshold, maximumResolution, initialResolution);
		}
		else
		{
			bool add = false;
			float maxValue = 0;
			for (int b = 0; b < 2; b++)
			{
				if ((*curchild)->neighborDifference(b, maxValue) > bandLevel)
				{
					add = true;
					break;
				}
			}

			//cout<<(*curchild)->leo_thr<<"\n";
			if (add && (*curchild)->leo_thr>0.0f)
			{

				targetx = ((*curchild)->x1 + (*curchild)->x2) / 2.0f;
				targety = ((*curchild)->y1 + (*curchild)->y2) / 2.0f;
				Connection2D* p;
				weight = scale((*curchild)->activationLevel, -5.0, 5.0);
				if (outgoing)
					p = new Connection2D((*curchild)->fixedx, (*curchild)->fixedy, targetx, targety, weight);
				else
					p = new Connection2D(targetx, targety, (*curchild)->fixedx, (*curchild)->fixedy, weight);
				list.push_back(p);
			}

		}

	}
}

void ESQuadTree::divide(int resolution)
{
	std::vector<ESQuadTree*> rectList;
	std::vector<ESQuadTree*> rectToAdd;

	//	List<Rect> rectList = new List<Rect>();
	//	List<Rect> rectToAdd = new List<Rect>();

	float midx, midy;
	rectList.push_back(this);
	//  minDistance = 0.25;
	int squares = 0;
	float x1, y1, x2, y2;

	std::vector<ESQuadTree*>::iterator it;

	int currentLevel = level;

	while (true)    //TODO replace. Do not calculate distance for every rectangle
	{
		//  (squares<(resolution*resolution))//
		rectToAdd.clear();

		if (currentLevel < resolution) //(Math.Abs(rectList[0].x1 - rectList[0].x2) > 0.4f))// || (sampleTreshold != -1.0f && dif > sampleTreshold))
		{
			for (it=rectList.begin();it!=rectList.end();it++)
			{
				//	cout<<(*it)->level<<"\n";
				midx = ( (*it)->x1 + (*it)->x2) / 2.0f;
				midy = ( (*it)->y1 + (*it)->y2) / 2.0f;

				rectToAdd.push_back(new ESQuadTree(fixedx, fixedy, outgoing, (*it)->x1,  (*it)->y1, midx, midy,  (*it)->cppn, (*it)));
				rectToAdd.push_back(new ESQuadTree(fixedx, fixedy, outgoing, (*it)->x1,  (*it)->y2, midx, midy,  (*it)->cppn, (*it)));
				rectToAdd.push_back(new ESQuadTree(fixedx, fixedy, outgoing, (*it)->x2,  (*it)->y1, midx, midy,  (*it)->cppn, (*it)));
				rectToAdd.push_back(new ESQuadTree(fixedx, fixedy, outgoing, (*it)->x2,  (*it)->y2, midx, midy,  (*it)->cppn, (*it)));
				//rect.visited = true;
			}

			rectList.clear();
			for (it=rectToAdd.begin();it!=rectToAdd.end();it++)
			{
				x1 = ((*it)->x2 + (*it)->x1) / 2.0f;
				y1 = ((*it)->y2 + (*it)->y1) / 2.0f;

				if (outgoing) //outgoing connectivity pattern
				{
					(*it)->activationLevel = queryCPPN(fixedx, fixedy, x1, y1,(*it)->leo_thr);
				}
				else	//incomming connectivity pattern
				{
					(*it)->activationLevel = queryCPPN(x1, y1, fixedx, fixedy,(*it)->leo_thr);
				}
				rectList.push_back( (*it) );
			}
			currentLevel += 1;
		}
		else
		{
			//	cout<<rectList.size()<<"\n";
			break;
		}
	}
}

float ESQuadTree::queryCPPN(float x1, float y1, float x2, float y2, float& thr)
{
	coordinates[0]=x1;
	coordinates[1]=y1;
	coordinates[2]=0.0f;
	coordinates[3]=x2;
	coordinates[4]=y2;
	coordinates[5]=0.0f;
	coordinates[6]=1.0f;

	cppn->flush();

	cppn->load_sensors(coordinates);
	//cppn->RecursiveActivation(); //TODO implement!!!

	for(int i = 0; i < 5; i++)
		cppn->activate();

	//genome.MultipleSteps(NETWORK_ITERATIONS);  //TODO query CPPN based on depth

	//thrs = 0.0f;
	thr = (float)cppn->outputs[3]->activation; //TODO 3?
	return (float)cppn->outputs[0]->activation;  //use weight
}

void ESQuadTree::getPoints(std::vector<float> &l)
{
	vector<ESQuadTree*>::iterator it;

	if (childs.size() > 0.0f)
	{
		for (it=childs.begin();it!=childs.end(); it++)
		{
			(*it)->getPoints(l);
		}
	}
	else
	{
		l.push_back(activationLevel);
	}
}

float ESQuadTree::maxHighestDif(int dimension, bool direction, float &maxValue)
{

	float v;

	//if (float.IsNaN(activationLevel))
	//{
	//	//activationLevel = queryCPPN((x1 + x2) / 2.0f, (y1 + y2) / 2.0f, (z1 + z2) / 2.0f, (w1 + w2) / 2.0f);
	//	Console.WriteLine("ERROR :" + x1 + " " + x2 + " " + y1 + " " + y2);

	//}

	float newx1 = x1, newx2 = x2, newy1 = y1, newy2 = y2;
	switch (dimension)
	{
	case 0:
		nextPos(newx1, newx2, x1, x2, direction);
		break;
	case 1:
		nextPos(newy1, newy2, y1, y2, direction);
		break;
	}

	float tx, ty, tz, tw;
	tx = (newx2 + newx1) / 2.0f;
	ty = (newy2 + newy1) / 2.0f; //(newx2 + newx1) / 2.0f; //(newy2 + newy1) / 2.0f;  //BUG was  
	if ((tx < -1.0) || (tx > 1.0) || (ty < -1.0) || (ty > 1.0))
	{
		//Outside of bounds
		return 0.0f;
	}
	float tmp;
	//v = queryCPPN(fixedx, fixedy, tx, ty);
	if (outgoing)
	{
		v = queryCPPN(fixedx, fixedy, tx, ty, tmp);
	}
	else
	{
		v = queryCPPN(tx, ty, fixedx, fixedy, tmp);
	}

	if (abs(v) > abs(maxValue)) maxValue = v;

	return abs(v - activationLevel);
}

void ESQuadTree::nextPos(float &new1, float &new2, float &v1, float &v2, bool direction)
{
	if (!direction)
	{
		new1 = v1 - abs(v1 - v2);
		new2 = v1;
	}
	else
	{
		new2 = v2 + abs(v1 - v2);
		new1 = v2;
	}
}

float ESQuadTree::neighborDifference(int dimension, float &maxValue)
{
	return min(maxHighestDif(dimension, true, maxValue), maxHighestDif(dimension, false, maxValue));        
}

float ESQuadTree::variance()
{
	float med;
	if (childs.size() == 0)
	{
		med = activationLevel;
		return 0.0f;
	}

	std::vector<float> l;
	getPoints(l);

	float m = 0.0f, v = 0.0f;
	std::vector<float>::iterator it;

	for(it=l.begin();it!=l.end();it++)
	{
		m += (*it);
	}
	m /= l.size();
	med = m;
	for(it=l.begin();it!=l.end();it++)
	{
		v += (float)pow( (*it) - m, 2);
	}
	v /= l.size();
	return v;
}

//bool lookup(const hash_set<const NNode*, nodehasher>& Set, NNode* node)
//{
//	hash_set<const NNode*, nodehasher>::const_iterator it = Set.find(node);
//	if (it!=Set.end()) node = (*it);
//	return (it==Set.end());
//}

void EvolvableSubstrate::find(std::vector<NNode*> &list, NNode* &node, bool& newNode)
{
	vector<NNode*>::iterator it;
	for (int i=0; i<list.size(); i++)
	{
		if (list[i]->x_sub_pos==node->x_sub_pos && list[i]->y_sub_pos==node->y_sub_pos)
		{
			delete node;
			node = list[i];
			newNode = false;
			return;

		}
	}
	newNode = true;

}

Network* EvolvableSubstrate::generateSubstrate(std::vector<NNode*> inputs, std::vector<NNode*> outputs, Network *cppn,
	float initialRes, float varianceThreshold, float bandingThreshold, float ESIterations, float divisionThreshold, int maximumRes, 
	float minX, float minY, float maxX, float maxY)// std::vector<Link*> &connections, std::vector<NNode*> &hiddenNodes)
{
	int sourceIndex, targetIndex, neuronCount = 0;
	float connectionCoordinates[4];

	//Keep a list of all discovered hidden nodes so we don't add them again
	//	std::vector<Point2D*> tabuList;
	//	std::vector<Point2D*> hiddenPos;

	std::vector<NNode*>::iterator curnode;

	std::vector<Connection2D*> connections2D;

	int sensorNeuronCount = inputs.size()+outputs.size();
	//TODO clear hidden neurons if there are any beforehand
	int hiddenNeuronCount=0;

	//	hash_set<const NNode*, nodehasher> hiddenNodes;
	vector<NNode*> hiddenNodes;

	bool newNode;

	//Discover hidden nodes that connect directly to the input 
	for(curnode=inputs.begin();curnode!=inputs.end(); curnode++)
	{
		ESQuadTree tree( (*curnode)->x_sub_pos, (*curnode)->y_sub_pos, true, minX, minY, maxX, maxY, cppn);// = new ESQuadTree(input.X, input.Y, true, null, 0, minX, minY, maxX, maxY, genome);

		tree.divide(initialRes);

		tree.getConnections(connections2D, varianceThreshold, bandingThreshold, divisionThreshold, maximumRes, initialRes);

		std::vector<Connection2D*>::iterator curConnection;
		for(curConnection=connections2D.begin();curConnection<connections2D.end();curConnection++)
		{
			NNode *targetNode = new NNode(NEURON, hiddenNeuronCount+sensorNeuronCount, HIDDEN);
			targetNode->x_sub_pos=(*curConnection)->endx;
			targetNode->y_sub_pos=(*curConnection)->endy;

			int idbefore = targetNode->node_id;

			//If the hidden neuron already exsists it is replaced here by the excisting one otherwise add it
			find(hiddenNodes, targetNode, newNode);

			//If it is new increase the hidden node count
			if (newNode)
			{
				hiddenNeuronCount++;
				hiddenNodes.push_back(targetNode);
			}
			else
			{
				if (idbefore==targetNode->node_id)
					cout<<"This shouldn't happen "<<targetNode->node_id;
			}
			//	std::cout<<targetNode->x_sub_pos<<" "<<targetNode->y_sub_pos<<" "<<(*curConnection)->weight<<"\n";

			Link* newlink=new Link( (*curConnection)->weight,(*curnode),targetNode, false);//,curlink->is_recurrent);
			//TODO check if recursion works

			//	cout<<(*curnode)->x_sub_pos<<" "<<(*curnode)->y_sub_pos<<" "<<targetNode->x_sub_pos<<"\t"<<targetNode->y_sub_pos<<"\n";

			(targetNode->incoming).push_back(newlink);
			(*curnode)->outgoing.push_back(newlink);

		}

		for(curConnection=connections2D.begin();curConnection<connections2D.end();curConnection++)
		{
			delete (*curConnection);
		}
		connections2D.clear();
	}
	//Discover hidden to hidden nodes
	//Discover hidden nodes that connect directly to the input 
	vector<NNode*> temp;
	temp.insert(temp.end(), hiddenNodes.begin(), hiddenNodes.end());

	for(curnode=hiddenNodes.begin();curnode!=hiddenNodes.end(); curnode++)
	{
		ESQuadTree tree( (*curnode)->x_sub_pos, (*curnode)->y_sub_pos, true, minX, minY, maxX, maxY, cppn);// = new ESQuadTree(input.X, input.Y, true, null, 0, minX, minY, maxX, maxY, genome);

		tree.divide(initialRes);

		tree.getConnections(connections2D, varianceThreshold, bandingThreshold, divisionThreshold, maximumRes, initialRes);

		std::vector<Connection2D*>::iterator curConnection;
		for(curConnection=connections2D.begin();curConnection<connections2D.end();curConnection++)
		{
			NNode *targetNode = new NNode(NEURON, hiddenNeuronCount+sensorNeuronCount, HIDDEN);
			targetNode->x_sub_pos=(*curConnection)->endx;
			targetNode->y_sub_pos=(*curConnection)->endy;

			//If the hidden neuron already exsists it is replaced here by the excisting one otherwise add it
			find(temp, targetNode, newNode);

			//If it is new increase the hidden node count
			if (newNode) 
			{
				temp.push_back(targetNode);
				hiddenNeuronCount++;
			}
			Link* newlink=new Link( (*curConnection)->weight,(*curnode),targetNode, false);//,curlink->is_recurrent);
			//TODO check if recursion works

			(targetNode->incoming).push_back(newlink);
			(*curnode)->outgoing.push_back(newlink);

		}

		for(curConnection=connections2D.begin();curConnection<connections2D.end();curConnection++)
		{
			delete (*curConnection);
		}
		connections2D.clear();
	}
	hiddenNodes.clear(); //FIXME probably not the smartest way
	hiddenNodes.insert(hiddenNodes.end(), temp.begin(), temp.end());

	//Hidden to Hidden
	//for (int step = 0; step < ESIterations; step++)
	//{
	//    pointToAdd.Clear();
	//    foreach (PointF hiddenP in hiddenPos)
	//    {
	//        Rect startRec = new Rect(hiddenP.X, hiddenP.Y, true, null, 0, minX, minY, maxX, maxY, genome);
	//        startRec.createTree(initialResolution);
	//        startRec.addPoints(ref _connections, ref varianceThreshold, ref bandThreshold, ref divsionThreshold, ref maximumResolution);

	//        foreach (ExpressPoint p in _connections)
	//        {
	//            targetX = (p.x1 + p.x2) / 2.0f;
	//            targetY = (p.y1 + p.y2) / 2.0f;
	//            PointF newp = new PointF(targetX, targetY);
	//            if (!tabuList.Contains(newp))
	//            {
	//                pointToAdd.Add(newp);
	//                tabuList.Add(newp);
	//            }
	//        }
	//    }
	//    hiddenPos.Clear();
	//    if (pointToAdd.Count == 0) break;
	//    hiddenPos.AddRange(pointToAdd);
	//}


	//foreach (ExpressPoint t in _connections)
	//{
	//    connectionCoordinates[0] = t.fixedx;
	//    connectionCoordinates[1] = t.fixedy;
	//    connectionCoordinates[2] = (float)(t.x1 + t.x2) / 2.0f;
	//    connectionCoordinates[3] = (float)(t.y1 + t.y2) / 2.0f;

	//    if (float.IsNaN(t.activationLevel))
	//    {
	//        Console.WriteLine("Normally this shouldn't happen");
	//        return;
	//        //  
	//    }
	//    else
	//    {
	//        output = t.activationLevel;
	//    }

	//    PointF source = new PointF(connectionCoordinates[0], connectionCoordinates[1]);
	//    PointF target = new PointF(connectionCoordinates[2], connectionCoordinates[3]);

	//    sourceIndex = hiddenNeurons.IndexOf(source);            //TODO change. computationally expensive
	//    if (sourceIndex == -1) //!hiddenNeurons.Contains(source)) 
	//    {
	//        sourceIndex = hiddenNeurons.Count;
	//        hiddenNeurons.Add(source);
	//        neuronCount++;
	//    }

	//    targetIndex = hiddenNeurons.IndexOf(target);
	//    if (targetIndex == -1) //!hiddenNeurons.Contains(target)) hiddenNeurons.Add(target);
	//    {
	//        targetIndex = hiddenNeurons.Count;
	//        hiddenNeurons.Add(target);
	//        neuronCount++;
	//    }

	//    weight = (float)(((Math.Abs(output) - (connectionThreshold)) / (1 - connectionThreshold)) * HyperNEATParameters.weightRange * Math.Sign(output));
	//    //if (weight > 0.0) weight = 1.0f;
	//    //else weight = -1.0f;

	//    connections.Add(new ConnectionGene(counter++, (uint)(sourceIndex + inputCount + outputCount), (uint)(targetIndex + inputCount + outputCount), weight, ref connectionCoordinates, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f));
	//    //   }
	//}

	//_connections.Clear();


	//Discover nodes connecting to  the outputs
	for(curnode=outputs.begin();curnode!=outputs.end(); curnode++)
	{
		ESQuadTree tree( (*curnode)->x_sub_pos, (*curnode)->y_sub_pos, false, minX, minY, maxX, maxY, cppn);// = new ESQuadTree(input.X, input.Y, true, null, 0, minX, minY, maxX, maxY, genome);

		tree.divide(initialRes);

		tree.getConnections(connections2D, varianceThreshold, bandingThreshold, divisionThreshold, maximumRes, initialRes);

		std::vector<Connection2D*>::iterator curConnection;
		for(curConnection=connections2D.begin();curConnection<connections2D.end();curConnection++)
		{
			NNode *sourceNode = new NNode(NEURON, hiddenNeuronCount+sensorNeuronCount, HIDDEN);
			sourceNode->x_sub_pos=(*curConnection)->startx;
			sourceNode->y_sub_pos=(*curConnection)->starty;

			//If the hidden neuron already exist it is replaced by the excisting one
			find(hiddenNodes, sourceNode, newNode);

			//If it is new increase the hidden node count
			if (!newNode) 
			{
				//Only add if this hidden node already exist otherwise it will not have a connection to the inputs
				Link* newlink=new Link( (*curConnection)->weight, sourceNode, (*curnode), false);//,curlink->is_recurrent);
				//TODO check if recursion works

				(*curnode)->incoming.push_back(newlink);
				sourceNode->outgoing.push_back(newlink);
				//cout<<"old";
			}


		}

		for(curConnection=connections2D.begin();curConnection<connections2D.end();curConnection++)
		{
			delete (*curConnection);
		}
		connections2D.clear();
	}

	//Get bias and time constant for hidden and outputs nodes	
	double coordinates[7], ur_bias, ur_time_const;

	for (curnode=hiddenNodes.begin(); curnode!=hiddenNodes.end();curnode++)
	{
		if (((*curnode)->type)!=SENSOR) {
			coordinates[0]=(*curnode)->x_sub_pos;
			coordinates[1]= (*curnode)->y_sub_pos;
			coordinates[2]= 0.0f; // Node's position
			coordinates[3]= 0.0f;
			coordinates[4]= 0.0;
			coordinates[5]= 0.0;
			coordinates[6]= 1.0f; 

			cppn->load_sensors(coordinates);

			// Activate the CPPN max_depth times to propagate the output through the entire network
			//for(int i = 0; i < max_depth; i++)
			for(int i = 0; i < 5; i++)
				cppn->activate();

			(*curnode)->ftype = ACTIVATION_FUNCTION_SIGMOID;
			ur_bias = cppn->outputs[1]->activation;
			ur_time_const = cppn->outputs[2]->activation;

			// Scale bias to [-8.0, 8.0]
			//ur_bias = scale(ur_bias, -8.0, 8.0);
			ur_bias = scale(ur_bias, -3.0, 3.0);
			//cout<<ur_bias; //TODO Risi: check if this is corret

			// Scale time constant within [0.1, 5.0]
			//ur_time_const = scale(ur_time_const, 0.1, 5.0);
			ur_time_const = scale(ur_time_const, 0.1, 2.0);

			(*curnode)->bias = ur_bias;  //TODO do we need to round?
			(*curnode)->time_const = ur_time_const;

			cppn->flush();
		}
	}

	vector<NNode*> allnodes;
	allnodes.insert(allnodes.end(),inputs.begin(),inputs.end());
	allnodes.insert(allnodes.end(),outputs.begin(),outputs.end());
	allnodes.insert(allnodes.end(),hiddenNodes.begin(),hiddenNodes.end());

//	pruneConnections(allnodes);
	Network* net = new Network(inputs,outputs,allnodes, 0); //TODO is the ID important
	return net;
}

void EvolvableSubstrate::pruneConnections(std::vector<NNode*> &list)
{
	bool danglingConnection = true;
	std::vector<NNode*>::iterator curnode;

	std::vector<Link*>::iterator curConnection;

	for(curnode=list.begin();curnode!=list.end(); curnode++)
	{
		//remove if no incomming connections
		if ( (*curnode)->incoming.size()==0 && (*curnode)->gen_node_label != INPUT)
		{
			for (curConnection=(*curnode)->outgoing.begin(); curConnection!= (*curnode)->outgoing.end(); curConnection++)
			{
				(*curConnection)->enabled = false;//in_node->outgoing.erase(curConnection);

			//	(*curConnection)->out_node->incoming.erase(curConnection);
			//	delete (*curConnection);
				//delete
			}
		//	delete (*curnode);
		}

		while (danglingConnection)
		{
			bool f=false;
			if ( (*curnode)->outgoing.size()==0  && (*curnode)->gen_node_label != OUTPUT)
			{
				for (curConnection=(*curnode)->incoming.begin(); curConnection!= (*curnode)->incoming.end(); curConnection++)
				{
					//(*curConnection)->enabled = false;//in_node->outgoing.erase(curConnection);
					f=true;
					(*curnode)->incoming.erase(curConnection);
					(*curConnection)->out_node->outgoing.erase(curConnection);
				//	delete (*curConnection);
					//delete
				}
			}
			if (!f) danglingConnection = false;
		}

	}
}