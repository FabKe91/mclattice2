#ifndef LIPIDPROPERTIES_H
#define LIPIDPROPERTIES_H
#include <vector>
#include <memory>
#include "inputfile.h"
#include <cmath>
#include <boost/multi_array.hpp>


class LipidProperties
{
public:
    LipidProperties();
    void readParas(std::shared_ptr<InputFile>);
    double** neighbourFunction;
    double** entropyFunction;
    double*** enthalpyFunction;
    
private:
    std::shared_ptr<InputFile> inputfile;
    
    
    double polynom(std::vector<double>&, double);
    double sigmoid(std::vector<double>&, double);


};

#endif // LIPIDPROPERTIES_H
