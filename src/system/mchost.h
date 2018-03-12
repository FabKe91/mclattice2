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
    int steps=0;
    int imageRate;
    int notAcceptedSwaps=0;
    int notAcceptedFlucs=0;
    int loopCounter=0;
    std::vector<int> IDs;



    std::shared_ptr<LipidProperties> lipidproperties;
    std::shared_ptr<InputFile> inputfile;
    std::unique_ptr<DataFile> datafile;
    
    void doSystemloop();

    double calcSwapEnthalpy();
    double calcHostFreeEnerg();
   
    int lastSwappedIDs[2]={0,0};    
    void setHost(int x, int y);
    void setHost(int ID);
    void setRNDHost();
    void setPartner();
    
    int rdnPartnerNumber=0;

    
    inline bool acceptance(const double Enthalpy1, const double Enthalpy2);

public:
    void run();   

    void setup(std::string);
    void setupForRestart(std::string);
    

};

#endif // MCHOST_H
