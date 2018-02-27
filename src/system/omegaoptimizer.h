#ifndef OMEGAOPTIMIZER_H
#define OMEGAOPTIMIZER_H
#include "lipidsystem.h"
#include <iostream>
#include <map>
#include <fstream>
#include "datafile.h"
#include "inputfile.h"
#include <chrono>
#include <ctime>
#include<iomanip>
#include<cmath>

class OmegaOptimizer
{
private:
    Lipidsystem lipidsystem;
    double T=1;
    double kB=1;
    double kBT=1;
    int steps=0;
    int imageRate;
    int notAcceptedSwaps=0;
    int notAcceptedFlucs=0;



    std::shared_ptr<LipidProperties> lipidproperties;


    void updateOmega(std::vector<int>);
    void runUntilEquilibrium();
    void doSystemloop();
    void calcCurrentOrderDestr(int);
    std::vector<double> MDOrderDestr;
    std::vector<double> currentOrderDestr;

    
public:
    void run();
    void optimizeOmega();

    void setupOptimization();
    
    bool acceptance(const double Enthalpy1, const double Enthalpy2);

};

#endif // OMEGAOPTIMIZER_H
