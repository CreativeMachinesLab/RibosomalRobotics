#ifndef POPUL123
#define POPUL123
#include <vector>
#include <algorithm>
#include "individual.h"
#include "helpers.h"

using namespace std;


namespace SGA
{
class Population {
    public:
       
       bool dirty;
	   Individual* prototype;
       Population()
       {
            prototype=NULL;
            dirty=false;
       }
       
       vector<Population*> subpopulations;
       vector<Individual*> individuals;
       virtual void generation_callback()
	   {
	   }
       virtual void insert(Individual* ind)
       {
        dirty=true;
        ind->incref();
        individuals.push_back(ind);
       }
       
       virtual void adjust_fitness()
       {
         for(int x=0;x<subpopulations.size();x++)
         {
            subpopulations[x]->adjust_fitness();
         } 
       }
       
       virtual void selfsort()
       {
            sort(individuals.begin(),individuals.end(),&fitcmp);
            dirty=false;
       }
       
       virtual void selfsort_raw()
       {
            sort(individuals.begin(),individuals.end(),&rawcmp);
            dirty=true;
       }
       
       virtual Individual* reproduce()
       {
            if (dirty)
                selfsort();
            int max=individuals.size();
            int min=(int)(0.8*max);         
            Individual* newind = individuals[rand_int(min,max)]->clone();
            newind->mutate();
            return newind;
       }
       
       virtual void remove_ind(Individual* ind)
       {
            vector<Individual*>::iterator it;
            for(it=individuals.begin();it!=individuals.end();it++)
            {
                if(*it==ind)
                {
                    individuals.erase(it);
                    ind->decref();
                    break;
                }
            }
       }
       
       virtual void remove_worst()
       {
            if (dirty)
                selfsort();
                
            Individual* worst=individuals.front();
            individuals.erase(individuals.begin());
            worst->decref(); 
            
            for(int x=0;x<subpopulations.size();x++)
                subpopulations[x]->remove_ind(worst);         
       }
       
       virtual void clear()
       {
            for(int x=0;x<individuals.size();x++)
                individuals[x]->decref();
                
            individuals.clear();
            
            for(int x=0;x<subpopulations.size();x++)
                subpopulations[x]->clear();
       }
};

class SpeciatedPopulation: public Population {
    
    public:
        int target_species;
        double compat_threshold;
        
    SpeciatedPopulation()
    {
        target_species = 5;
        compat_threshold = 0.2;
    }
    
	virtual void generation_callback()
	{
		if (subpopulations.size()<target_species)
		{
			compat_threshold*=0.9;
		}
		else if (subpopulations.size()>target_species)
		{
			compat_threshold*=1.1;
		}
	}
    virtual void adjust_fitness()
    {
        int size=subpopulations.size();
        for(int x=0;x<size;x++)
        {
            int subsize=subpopulations[x]->individuals.size();
            
            //fitness sharing
            for(int y=0;y<subsize;y++)
                subpopulations[x]->individuals[y]->fitness->fitness /= subsize;
        }
        Population::adjust_fitness();
    }
    
    
    
    
    virtual void add_new_species(Individual* proto)
    {
        Population* newpop = new Population();
        newpop->prototype=proto;
		proto->incref();
        newpop->insert(proto);
	
		subpopulations.push_back(newpop);
    }
    
    virtual void insert(Individual* ind)
    {
        bool match=false;
        Population* matchpop;
        for(int x=0;x<subpopulations.size();x++)
        {
            if(ind->distance(subpopulations[x]->prototype)<compat_threshold)
            {
                //break on first match
                match=true;
                matchpop=subpopulations[x];
                break;
            }    
        }
        
        if(!match)
        {
			cout << "add new species" << endl;
            add_new_species(ind);
        }
        else
        {
			//cout << "insert" << endl;
            matchpop->insert(ind);
        }
        
        //call superclass
        Population::insert(ind);
    }
    
    virtual void remove_worst()
    {
        Population::remove_worst();
        vector<Population*>::iterator it;
        for(it=subpopulations.begin();it!=subpopulations.end();it++)
        {
            //delete species if no members
            if((*it)->individuals.size()==0)
            {
				cout << "delete species" << endl;
				(*it)->prototype->decref();
				delete (*it);
				subpopulations.erase(it);
				it=subpopulations.begin();
			}
        }
    }
            
};
}
#endif

