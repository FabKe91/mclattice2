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
    void readParas(std::shared_ptr<LipidProperties>);
    void setup();
   
    
    unsigned int getHeight(){ return height; };
    unsigned int getWidth(){ return width; };
    
    
    const boost::multi_array<int,2> getOrderParas();
    const boost::multi_array<int,2> getTypes();
    const boost::multi_array<int,2> getIDs();
    int getMeanOrder();
    std::vector<int> getOrderDestr();

    double calcSwapEnthalpy();
    double calcHostFreeEnerg();

    void findPair();
    void swap();
    
    void setHost(int x, int y);
    void setHost(int ID);
    void setRNDHost();
    void setPartner();
   
    
    
    void fluctuate();
    void fluctuateBack();
    
    std::shared_ptr<InputFile> inputfile;
private:
    std::shared_ptr<LipidProperties> lipidproperties;
    unsigned int height;
    unsigned int width;
    int lastSwappedIDs[2]={0,0};

    int oldOrder=0;
    std::vector<Lipid> lipids {};
    

    int** map;
    int rdnPartnerNumber=0;
    
    void printMap();//only for debugging

    
    double calcPairEnthalpy(unsigned int ID1,unsigned int ID2);


};

#endif // LIPIDSYSTEM_H
