#ifndef INPUTFILE_H
#define INPUTFILE_H
#include<fstream>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include "typeproperties.h"
class InputFile
{
    
public:
    double T, kBT;
    int width,height;
    int nType=0; //counting types 
    std::map<std::string,double> paras; //general parameter map
    std::map<std::string,int> typeMap; //mapping typenames to typeIDs
    std::vector<std::string> outs; //strings to print (see datafile)
    std::vector<TypeProperties> types;
    std::vector<double> concentrations;
    
    std::vector<std::vector<double>> entropyPara;
    std::vector<std::vector<double>> selfEnergiePara;
    std::vector<std::vector<std::vector<double>>> lipidCholEnergiePara;
    std::vector<std::vector<std::vector<std::vector<double>>>> enthalpyPara;
    
    
    std::vector<std::vector<double>>  cholLipidNeighPara;
    std::vector<std::vector<std::vector<double>>> LipidLipidNeighPara;
    std::vector<double> CholCholEnergiePara;
    
    InputFile(std::string);
};

#endif // INPUTFILE_H
