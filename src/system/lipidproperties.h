#ifndef LIPIDPROPERTIES_H
#define LIPIDPROPERTIES_H
#include <vector>
#include <memory>
#include "inputfile.h"
#include <cmath>
#include <boost/multi_array.hpp>
#include "../lib/enhance.hpp"
#include "inputfile.h"




class LipidProperties
{
public:
    LipidProperties();
    void readParas();
    double** neighbourFunction;
    double** entropyFunction;
    double** selfEnergieFunction;
    double*** enthalpyFunction;
    
    
    void updateKBT();
    
    double polynom(std::vector<double>&, double);
    double sigmoid(std::vector<double>&, double);

private:
    
    


};

#endif // LIPIDPROPERTIES_H
