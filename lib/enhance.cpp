#include "enhance.hpp"


namespace enhance 
{
    unsigned int    seed;
    std::mt19937_64 rand_engine;


//     // random double from [a,b)
//     double random_double(double a, double b)
//     {
//         std::uniform_real_distribution<double> distribution(a,b);
//         return distribution(rand_engine);
//     }
// 
//     // random int from [a,b]
//     unsigned int random_uns_int(int a, int b)
//     {
//         std::uniform_int_distribution<int> intdistribution(a,b);
//         return intdistribution(rand_engine);
//     }
// 
//     
//     // random int from [a,b]
//     int random_int(int a, int b)
//     {
//         std::uniform_int_distribution<int> intdistribution(a,b);
//         return intdistribution(rand_engine);
//     }
// 
// 
//     float fastExp(float x)
//     {
//         x = 1.0 + x / 256.0;
//         x *= x; x *= x; x *= x; x *= x;
//         x *= x; x *= x; x *= x; x *= x;
//         return x;        
//     }

    double polynom(std::vector<double>& coeff, double x)
    {
        double y=0;
        for(int n=0; n<coeff.size();n++)
        {
            y+=std::pow(x,n)*coeff[n];
        }
        return y;
    }

    double sigmoid(std::vector<double>& coeff, double T)
    {
        if(coeff.size()!=4) throw std::invalid_argument("need 4 parameters for sigmoid");
        return coeff[0] + (coeff[1] / ( 1+std::exp(-coeff[2]*(T-coeff[3])) ) );  
    }

    double logistic(std::vector<double>& coeff, double S)
    {
        if(coeff.size()!=4) throw std::invalid_argument("need 4 parameters for logistic");
        return coeff[0] + (coeff[1] / ( 1+std::exp(-coeff[2]*(S-coeff[3])) ) );  
    }
    
    double lennard(std::vector<double>& coeff, double S)
    {
        if(coeff.size()!=4) throw std::invalid_argument("need 4 parameters for lennard");
        double Q = std::pow( ( (coeff[0] - coeff[1]) / (S - coeff[1]) ), coeff[2]);
        return -coeff[3] * Q * (Q-2);
    }
}
