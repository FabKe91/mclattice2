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

    double random_double(double, double);
    unsigned int    random_uns_int(int, int);
    int    random_int(int, int);

    float fastExp(float x);
    

}


namespace enh = enhance;
