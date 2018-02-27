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
    static int nType; //counting types 
    static std::map<std::string,double> paras; //general parameter map
    static std::map<std::string,int> typeMap; //mapping typenames to typeIDs
    static std::vector<std::string> outs; //strings to print (see datafile)
    static std::vector<TypeProperties> types;
    static std::vector<double> concentrations;
    
    static std::vector<std::vector<double>> neighbourPara;
    static std::vector<std::vector<double>> entropyPara;
    static std::vector<std::vector<double>> selfEnergiePara;
    static std::vector<std::vector<std::vector<double>>> enthalpyPara;
    
    static void loadFile(std::string);
};

#endif // INPUTFILE_H
