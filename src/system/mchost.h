#ifndef MCHOST_H
#define MCHOST_H
#include "lipidsystem.h"
#include <iostream>
#include <map>
#include <fstream>
#include "../IO/datafile.h"
#include "../IO/inputfile.h"
#include "cholesterinsystem.h"
#include <chrono>
#include <ctime>

class MCHost
{
protected:
    Lipidsystem lipidsystem;
    CholesterinSystem cholesterinsystem;
    int steps=0;
    int imageRate;
    int notAcceptedSwaps=0;
    int notAcceptedFlucs=0;
    int loopCounter=0;
    int CholSwaps=0;
    int notAcceptedCholSwaps=0;
    std::vector<int> IDs;
    std::vector<int> cholIDs;

    


    std::shared_ptr<LipidProperties> lipidproperties;
    std::shared_ptr<InputFile> inputfile;
    std::unique_ptr<DataFile> datafile;
    
    void doSystemloop();

    double calcSwapEnthalpy();
//     double calcEnthalpyAfterSwap();
    double calcHostFreeEnerg();
   
    int hostID=0;
    std::array<int,4> hostNeighbours;
    
    int partnerID=0;
    std::array<int,4> partnerNeighbours;
    
    
    void setHost(int x, int y);
    void setHost(int ID);
    void setRNDHost();
    void setPartner();
    
    int findLipidPairCholNeighbours(int ID1, int ID2);
    int findCholPairCholNeighbours(int ID1, int ID2);
    int findLipidCholPairCholNeighbours(int, int );
    
    
    inline bool acceptance(const double Enthalpy1, const double Enthalpy2);
    
    
    bool setRNDCholHost();
    bool setCholHost(int);
    bool setCholPartner();
    
    int cholHostID,cholPartnerID;
    
    double calcCholSwapEnerg();
    
    std::array<int,4> getLipidNeighOfChol(int);
    std::array<int,4> getLipidNeighOfLipid(int);
    std::array<int,4> getCholNeighOfChol(int);
    std::array<int,4> getCholNeighOfLipid(int);
    int getNumberCholNeighOfLipid(int ID);
    int getNumberCholNeighOfChol(int ID);

    
    bool IDinArrayLen4(int& ID, std::array<int,4>&);
   
    
    

public:
    void run();   

    void setup();
    void setupForRestart();
    

};

#endif // MCHOST_H
