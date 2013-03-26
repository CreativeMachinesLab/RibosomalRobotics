#ifndef INDIV123
#define INDIV123
#include <vector>
#include "genome.h"
#include "fitness.h"
using namespace std;


namespace SGA
{
class Individual {
    public:
    int ref_count;
    Individual(GenomePart* g, Fitness* f)
    {
        genome=g;
        fitness=f;
        ref_count=0;
    }
    void incref()
    {
        ref_count++;
    }
    void decref()
    {
        ref_count--;
        if (ref_count <= 0)
            delete this;
    }
    
    virtual Individual* clone()
    {
        GenomePart* g=genome->clone();
        Fitness* f=fitness->clone();
        Individual* newind=new Individual(g,f);
        return newind;
    }

    ~Individual()
    {
        delete genome;
        delete fitness;
    }
    GenomePart* genome;
    Fitness* fitness;
    virtual void mutate()
    {
        genome->mutate();
    }
    virtual double distance(Individual* ind)
    {
        return genome->distance(ind->genome);
    }

};


 bool rawcmp(const Individual *a , const Individual *b)
 {
   return  a->fitness->raw_fitness < b->fitness->raw_fitness;
 }
 
 bool fitcmp(const Individual *a , const Individual *b)
 {
    return a->fitness->fitness < b->fitness->fitness;
 }

typedef Individual* (*Generator)();
}
#endif

