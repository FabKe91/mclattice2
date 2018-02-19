#ifndef MCHOST_H
#define MCHOST_H
#include "lipidsystem.h"
#include <iostream>
#include <map>
#include <fstream>
#include "datafile.h"
#include "inputfile.h"
#include <chrono>
#include <ctime>

class MCHost
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



    std::shared_ptr<InputFile> inputfile;
    std::shared_ptr<LipidProperties> lipidproperties;
    std::unique_ptr<DataFile> datafile;
    
public:
    void run();
    void optimizeOmega();

    void setup(std::string);
    bool acceptance(const double Enthalpy1, const double Enthalpy2);

};

#endif // MCHOST_H
