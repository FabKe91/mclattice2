#ifndef LIPIDPROPERTIES_H
#define LIPIDPROPERTIES_H
#include <vector>
#include <memory>
#include "../IO/inputfile.h"
#include <cmath>
#include <boost/multi_array.hpp>
#include "../lib/enhance.hpp"




class LipidProperties
{
public:
    LipidProperties();
    void readParas(std::shared_ptr<InputFile>);
    double*** neighbourFunction;
    double** entropyFunction;
    double** selfEnergieFunction;
    double**** enthalpyFunction;
    double*** lipidCholEnergieFunction;
    double* cholCholEnergie;
    double** cholLipidNeigh; //N^CL(N^CC)
    
    
    void updateKBT();
    
    double polynom(std::vector<double>&, double);
    double sigmoid(std::vector<double>&, double);

private:
    std::shared_ptr<InputFile> inputfile;

    
    


};

#endif // LIPIDPROPERTIES_H
