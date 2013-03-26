#ifndef FITNESS123
#define FITNESS123

namespace SGA
{
class Fitness {
    public:
    double fitness;
    double raw_fitness;
    virtual Fitness* clone()
    {
        Fitness* new_fit=new Fitness();
        new_fit->fitness=fitness;
        new_fit->raw_fitness=raw_fitness;
        return new_fit;
    }
    bool virtual operator<(const Fitness &a)  
    {
        return fitness < a.fitness;
    }    
};
}
#endif


