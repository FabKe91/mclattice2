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

    /**
     * NOTE: Primarily used for neighbor rescaling of DPPC (function of temperature)
     * NOTE: This function is similar to sigmoid, but introduced to discern from sigmoid used for neighbor scaling function
     * logistic function of form:
     *    f(S) = Nmax + ( deltaN / (1+exp(-k(T-transitionT))) )
     * @param[in] coeff list of coeffs: <Nmax> <deltaN> <k> <transitionT>
     * @param[in] T function variable
     */
    double sigmoid(std::vector<double>& coeff, double T);

    /**
     * @param[in] coeff list of coefficients for parametrizing polynomial of degree coeff.size-1 (ascending order)
     * @param[in] x function variable
     */
    double polynom(std::vector<double>& coeff, double x);

    /**
     * NOTE: Primarily used for PL-PL interactions (function of order)
     * NOTE: This function is similar to sigmoid, but introduced to discern from sigmoid used for neighbor scaling function
     * logistic function of form:
     *    f(S) = Emax + ( Egain / (1+exp(-k(S-S0))) )
     * @param[in] coeff list of coeffs: <Emax> <Egain> <k> <S0>
     * @param[in] S function variable
     */
    double logistic(std::vector<double>& coeff, double x);

    /**
     * NOTE: Primarily used for PL-CHOL interaction
     * lennard jones function of form:
     *      Q = [(S0 - Sshift) / (S - Sshift)]**k
     *      f(S) = -E0 * Q * (Q - 2)
     * @param[in] coeff list of coeffs: <S0> <Sshift> <k> <E0>
     * @param[in] S function variable
     */
    double lennard(std::vector<double>& coeff, double x);




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
