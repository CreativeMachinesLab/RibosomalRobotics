#ifndef GA123
#define GA123
#include "population.h"
#include "fitness.h"
#include "genome.h"
#include "individual.h"
#include "population.h"
#include "domain.h"
#include "helpers.h"

#include <vector>
using namespace std;


namespace SGA
{

class GA {
    public:
        int epoch_counter;
        Population* population;
        Domain* domain;
        Generator gen;
        GA() {
            epoch_counter=0;
             }
        GA(Population* pop,Domain* dom,Generator g,int size)
        {
            population=pop;
            domain=dom;
            gen=g;
            for(int i=0;i<size;i++)
            {
                pop->insert(g());
            }
          epoch_counter=0;
          evaluate();
        }
		double best_rawfit()
		{
			double best=0.0;
			for (int x=0;x<population->individuals.size();x++)
			{
				double next=population->individuals[x]->fitness->raw_fitness;
				if (next > best)
					best=next;
			}
			return best;
		}

        double best_fit()
        {
            population->selfsort();
            int size=population->individuals.size();

            return population->individuals[size-1]->fitness->fitness;
        }
        
        virtual void evaluate()
        {
            int size=population->individuals.size();
            for(int x=0;x<size;x++)
                domain->evaluate(population->individuals[x]);
            population->adjust_fitness();
            population->selfsort();
        }
        
        virtual void epoch()
        {
            vector<Individual*> temp_pop;
            int orig_size=population->individuals.size();
			
			//species elitism
			for(int x=0;x<population->subpopulations.size();x++)
				if(population->subpopulations[x]->individuals.size()>5)
					temp_pop.push_back(population->subpopulations[x]->individuals[0]->clone());

            for(int x=0;x<orig_size-temp_pop.size();x++)
                temp_pop.push_back(population->reproduce());

            population->clear();
            for(int x=0;x<orig_size;x++)
            {
                domain->evaluate(temp_pop[x]);
                population->insert(temp_pop[x]);
            }
            population->adjust_fitness();
            epoch_counter+=1;
			population->generation_callback();
        }
};

class SteadyGA: public GA {
    public:
        SteadyGA(Population* pop, Domain* dom,Generator g,int size):GA(pop,dom,g,size)
        {

        }
        virtual void epoch()
        {
            int size=population->individuals.size();     
      
            for(int x=0;x<size;x++)
            {
                Individual* new_ind = population->reproduce();
                domain->evaluate(new_ind);
                population->insert(new_ind);
                population->adjust_fitness();
                population->remove_worst();
            }
            epoch_counter+=1;
			population->generation_callback();
        }
};
}
#endif

