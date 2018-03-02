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
    int steps=0;
    int imageRate;
    int notAcceptedSwaps=0;
    int notAcceptedFlucs=0;
    std::vector<int> IDs;
    int type;




    std::shared_ptr<LipidProperties> lipidproperties;
    std::shared_ptr<InputFile> inputfile;


    void updateOmega(std::vector<int>);
    void runUntilEquilibrium();
    void doSystemloop();
    void calcCurrentOrderDistr(int);
    std::vector<double> MDOrderDistr;
    std::vector<double> currentOrderDistr;

    
public:
    void run();
    void optimizeOmega();

    void setupOptimization(std::string,std::string);
    
    bool acceptance(const double Enthalpy1, const double Enthalpy2);

};

#endif // OMEGAOPTIMIZER_H
