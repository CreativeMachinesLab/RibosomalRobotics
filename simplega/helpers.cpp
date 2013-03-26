#include <cstdlib>
double rand_float()
{
    return rand() / (static_cast<double>(RAND_MAX));
}

double rand_float(float min, float max)
{
	return (max-min)*rand_float()+min;
}

int rand_int(int min,int max)
{
    return (rand() % (max-min)) + min;
}
