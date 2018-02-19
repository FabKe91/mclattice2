#include "enhance.hpp"


namespace enhance 
{
    unsigned int    seed;
    std::mt19937_64 rand_engine;


    // random double from [a,b)
    double random_double(double a, double b)
    {
        std::uniform_real_distribution<double> distribution(a,b);
        return distribution(rand_engine);
    }

    // random int from [a,b]
    unsigned int random_uns_int(int a, int b)
    {
        std::uniform_int_distribution<int> intdistribution(a,b);
        return intdistribution(rand_engine);
    }

    
    // random int from [a,b]
    int random_int(int a, int b)
    {
        std::uniform_int_distribution<int> intdistribution(a,b);
        return intdistribution(rand_engine);
    }


    float fastExp(float x)
    {
        x = 1.0 + x / 256.0;
        x *= x; x *= x; x *= x; x *= x;
        x *= x; x *= x; x *= x; x *= x;
        return x;        
    }
}
