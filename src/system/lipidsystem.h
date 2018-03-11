#ifndef LIPIDSYSTEM_H
#define LIPIDSYSTEM_H
#include "lipid.h"
#include <vector>
#include "../lib/enhance.hpp"
#include <iostream>
#include <fstream>
#include <map>
#include <boost/multi_array.hpp>
#include "inputfile.h"
#include "lipidproperties.h"
#include <cmath>

 
class Lipidsystem
{
public:
    Lipidsystem();
    void readParas(std::shared_ptr<LipidProperties>,std::shared_ptr<InputFile>);
    void setup();
   
    
   
    
    const boost::multi_array<int,2> getOrderParas();
    const boost::multi_array<int,2> getTypes();
    const boost::multi_array<int,2> getIDs();
    int getMeanOrder();
    std::vector<int> getOrderDistr();

    double calcSwapEnthalpy();
    double calcHostFreeEnerg();

    void swap();
    
    void setHost(int x, int y);
    void setHost(int ID);
    void setRNDHost();
    void setPartner();

    void setTypes(boost::multi_array<int,2>);
    void setOrder(boost::multi_array<int,2>);
    
    
    void fluctuate();
    void fluctuateBack();
    
private:
    std::shared_ptr<LipidProperties> lipidproperties;
    std::shared_ptr<InputFile> inputfile;

    unsigned int height;
    unsigned int width;
    int lastSwappedIDs[2]={0,0};

    int oldOrder=0;
    std::vector<Lipid> lipids {};
    

    int** map;
    int rdnPartnerNumber=0;
    
    void printMap();//only for debugging

    
    inline double calcPairEnthalpy(unsigned int ID1,unsigned int ID2);


};

#endif // LIPIDSYSTEM_H
