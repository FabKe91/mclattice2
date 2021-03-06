#pragma once

#include <ostream>
#include <cstdint>
#include <algorithm>
#include <random>
#include <iterator>
#include <type_traits>
#include <string>
#include <sys/stat.h>

namespace enhance
{

    extern unsigned int     seed;
    extern std::mt19937_64  rand_engine;

    inline double random_double(double, double);
    inline unsigned int    random_uns_int(int, int);
    inline int    random_int(int, int);

    inline float fastExp(float x);
    double sigmoid(std::vector<double>& coeff, double x);
    double polynom(std::vector<double>& coeff, double x);

    



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
namespace enh = enhance;
