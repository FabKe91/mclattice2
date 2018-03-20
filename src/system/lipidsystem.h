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

    void swap(int ID_0, int ID_1);
    

    void setTypes(boost::multi_array<int,2>);
    void setOrder(boost::multi_array<int,2>);
    
    
    void fluctuate(int ID);
    void fluctuateBack(int ID);

    inline double calcPairEnthalpy(int ID_1,int ID_2, int cholNeighs1,  int cholNeighs2, int pairCholNeighs);
    
    int** map;

    std::vector<Lipid> lipids {};

private:
    std::shared_ptr<LipidProperties> lipidproperties;
    std::shared_ptr<InputFile> inputfile;

    int oldOrder=0;
    


    
    void printMap();//only for debugging

    


};

double Lipidsystem::calcPairEnthalpy(int ID_1,int ID_2, int cholNeighs1,  int cholNeighs2, int pairCholNeighs)
{ 
    #ifndef NDEBUG
    std::cout<<"Lipidsystem::calcPairEnthalpy"<<std::endl;
    #endif
    
    
    int type1=lipids[ID_1].getType();
    int type2=lipids[ID_2].getType();
    int order1=lipids[ID_1].getOrder();
    int order2=lipids[ID_2].getOrder();
    

    if (type1>=type2) 
    {
        return lipidproperties->enthalpyFunction[type1][type2][pairCholNeighs][(int)((order1+order2)/2)]*(lipidproperties->neighbourFunction[type1][cholNeighs1]+lipidproperties->neighbourFunction[type2][cholNeighs2])/16;
    }
    else
    {
        return lipidproperties->enthalpyFunction[type2][type1][pairCholNeighs][(int)((order1+order2)/2)]*(lipidproperties->neighbourFunction[type1][cholNeighs1]+lipidproperties->neighbourFunction[type2][cholNeighs2])/16;;
    }

}

#endif // LIPIDSYSTEM_H
