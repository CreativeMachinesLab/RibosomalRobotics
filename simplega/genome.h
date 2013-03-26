#include <vector>
#include "helpers.h"
#ifndef GENOME123
#define GENOME123
using namespace std;

namespace SGA
{
class GenomePart {
    public:
    vector<GenomePart*> components;
    virtual GenomePart* clone()=0;
    virtual void mutate()=0;
    virtual double distance(GenomePart* p)=0;
    int type;
};

class IntegerGenome: public GenomePart {
    public:
        vector<int> vals;
        int max,min;
        IntegerGenome()
        {
        }
        IntegerGenome(int size,int min_ini, int max_ini)
        {
            max=max_ini;
            min=min_ini;
            for(int i=0;i<size;i++)
                vals.push_back(rand_int(min,max));
        }
        virtual void mutate()
        {
            for(int x=0;x<vals.size();x++)
            {
                if (rand_float()<0.2)
                {
                    vals[x]+=(rand_int(0,6)-3);
                    if (vals[x]>max)
                        vals[x]=max;
                    if (vals[x]<min)
                        vals[x]=min;
                }
            }
        }
        virtual double distance(GenomePart* p)
        {
            double accum=0.0;
            IntegerGenome* k = (IntegerGenome*)p;
            for(int x=0;x<vals.size();x++)
                accum+=abs(k->vals[x]-vals[x]);
            return accum;
        }
        GenomePart* clone()
        {
            IntegerGenome* the_clone= new IntegerGenome();
            the_clone->min=min;
            the_clone->max=max;
            for(int x=0;x<vals.size();x++)
                the_clone->vals.push_back(vals[x]);
            return (GenomePart*)the_clone;
        }
};

class FloatGenome: public GenomePart {
    public:
        vector<double> vals;
        int max,min;
		double mutsize;
        FloatGenome()
        {
        }
        FloatGenome(int size,double min_ini, double max_ini)
        {
            max=max_ini;
            min=min_ini;
			mutsize=(max-min)/10.0;
            for(int i=0;i<size;i++)
                vals.push_back(rand_float(min,max));
        }
        virtual void mutate()
        {
            for(int x=0;x<vals.size();x++)
            {
                if (rand_float()<0.2)
                {
                    vals[x]+=(rand_float(-mutsize,mutsize));
                    if (vals[x]>max)
                        vals[x]=max;
                    if (vals[x]<min)
                        vals[x]=min;
                }
            }
        }
        virtual double distance(GenomePart* p)
        {
            double accum=0.0;
            FloatGenome* k = (FloatGenome*)p;
            for(int x=0;x<vals.size();x++)
                accum+=abs(k->vals[x]-vals[x]);
            return accum;
        }
        GenomePart* clone()
        {
            FloatGenome* the_clone= new FloatGenome();
            the_clone->min=min;
            the_clone->max=max;
			the_clone->mutsize=mutsize;
            for(int x=0;x<vals.size();x++)
                the_clone->vals.push_back(vals[x]);
            return (GenomePart*)the_clone;
        }
};
}
#endif

