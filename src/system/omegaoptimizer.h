#ifndef OMEGAOPTIMIZER_H
#define OMEGAOPTIMIZER_H
#include "mchost.h"

class OmegaOptimizer : public MCHost
{
private:
    int type;



    void updateOmega(std::vector<int>);
    void runUntilEquilibrium();
    void calcCurrentOrderDistr(int);
    std::vector<double> MDOrderDistr;
    std::vector<double> currentOrderDistr;

public:
    void optimizeOmega();
    void setupOptimization(std::string);
    

};

#endif // OMEGAOPTIMIZER_H
