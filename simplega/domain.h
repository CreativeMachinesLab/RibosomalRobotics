#ifndef DOMAINNNN
#define DOMAINNNN
#include "genome.h"
#include <math.h>

namespace SGA
{
class Domain {
    public:
        virtual void evaluate(Individual* ind)=0;
};

class SimpleIntDomain:public Domain
{
public:
    virtual void evaluate(Individual* ind)
    {
        double accum=0.0;
        IntegerGenome* g=(IntegerGenome*)ind->genome;
        for(int x=0;x<g->vals.size();x++)
        {
            accum+=g->vals[x];
        }
        ind->fitness->fitness=accum;
        ind->fitness->raw_fitness=accum;
    }
};

class SimpleFloatDomain:public Domain
{
public:
    virtual void evaluate(Individual* ind)
    {
        double x=0.0,accum=0.0;
        FloatGenome* g=(FloatGenome*)ind->genome;
        x=g->vals[0];

		accum=exp(-2.0*log(2.0)*pow((x-0.1)/0.8,2.0))*pow(sin(3.14*5.0*x),6);
        
		ind->fitness->fitness=accum;
        ind->fitness->raw_fitness=accum;
    }
};
}
#endif

