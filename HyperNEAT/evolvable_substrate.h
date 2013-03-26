#ifndef _EVOLVABLE_SUBSTRATE_
#define _EVOLVABLE_SUBSTRATE_

#include "network.h"
#include "link.h"

namespace NEAT {
	
	struct Point2D
	{
		float x, y;
	};

	//A connection in space
	struct Connection2D
	{
		float startx, starty, endx, endy;
    //    float fixedx, fixedy;
        float weight;

        Connection2D(float _startx, float _starty, float _endx, float _endy, float _weight) :
		weight(_weight),startx(_startx),starty(_starty),endx(_endx),endy(_endy)
            {
            }
	};
		
	//The basic datastructure used by ES-HyperNEAT
	//...
	class ESQuadTree
	{
	
	public:
		ESQuadTree(float _fixedx, float _fixedy, bool _outgoing, float _x1, float _y1, float _x2, float _y2, Network* _cppn, ESQuadTree* _parent=NULL):
		  fixedx(_fixedx), fixedy(_fixedy),parent(_parent), outgoing(_outgoing), cppn(_cppn)
		{
               // activationLevel = float.NaN;
               // visited = false;
               // childs = new List<Rect>();
    
            if (parent != NULL)
            {
				parent->childs.push_back(this);
				level = parent->level+1;
				
            }
			else
			{
				level = 0;
			}


            if (_x1 < _x2)
            {
                x1 = _x1;
                x2 = _x2;
            }
            else
            {
                x2 = _x1;
                x1 = _x2;
            }

            if (_y1 < _y2)
            {
                y1 = _y1;
                y2 = _y2;
            }
            else
            {
                y2 = _y1;
                y1 = _y2;
            }
		}
		~ESQuadTree();

	    void divide(int resolution);
		void getConnections(std::vector<Connection2D*> &list, const float &treshold, const float &bandLevel, const float &divisionThreshold, const int &maximumResolution, const int &initialResolution);

	private:
		 double coordinates[7];
		 Network* cppn;
		 float x1, y1, x2, y2;
		 ESQuadTree *parent;
		 std::vector<ESQuadTree*> childs;
		 bool outgoing;
 		 float fixedx, fixedy, activationLevel; 
 		 int level; //in the quadtree
		 float leo_thr;
	
		 void nextPos(float &new1, float &new2, float &v1, float &v2, bool direction);
	     float neighborDifference(int dimension, float &maxValue);
		 float variance();
	     float maxHighestDif(int dimension, bool direction, float &maxValue);
         float queryCPPN(float x1, float y1, float x2, float y2, float& thr);
		 void getPoints(std::vector<float> &l);
	};

	class EvolvableSubstrate
	{

	public:
		//QuadTree* quadtree;

		EvolvableSubstrate();
		//Destructor 
		~EvolvableSubstrate() {} //TODO make sure to delete everything at the end

		//Generate the connections and hidden nodes
		Network* generateSubstrate(std::vector<NNode*> inputs, std::vector<NNode*> outputs, Network *cppn,
			float initialRes, float varianceThreshold, float bandingThreshold, float ESIterations, float divisionThreshold, int maximumRes, 
			float minX, float minY, float maxX, float maxY);

	private:
		void find(std::vector<NNode*> &list, NNode* &node, bool &newNode);
		void pruneConnections(std::vector<NNode*> &list);
	};

	//Helper class
	//struct eqnode
	//{
	//  bool operator()(const NNode* s1, const NNode* s2) const
	//  {
	//	  return (s1->x_sub_pos==s2->x_sub_pos) && (s1->y_sub_pos==s2->y_sub_pos);
	//  }
	//};

	//class nodehasher : public stdext::hash_compare <NNode>
	//{
	//	public:
	//	  size_t operator() (const NNode& p) const
	//	  // This is the hash function itself
	//	  {
	//		  return p.x_sub_pos * 100 + p.y_sub_pos;
	//	  }

	//	  bool operator() (const NNode& s1, const NNode& s2) const
	//	  // This is a comparison operator. It is needed by the hashed container
	//	  {
	//		 return (s1.x_sub_pos==s2.x_sub_pos) && (s1.y_sub_pos==s2.y_sub_pos);
	//	  }
	//};
}
#endif