#ifndef CHOLESTERINSYSTEM_H
#define CHOLESTERINSYSTEM_H

#include "../IO/inputfile.h"
#include "cholesterin.h"
#include "../lib/enhance.hpp"

#include <vector>
#include <algorithm>
#include <stdio.h>
#include <iostream>
#include <memory>
#include <boost/multi_array.hpp>


class CholesterinSystem
{
private:
    std::shared_ptr<InputFile> inputfile;


public:
    int** map;

    std::vector<Cholesterin> chols {};

    CholesterinSystem();
    void setup(std::shared_ptr<InputFile>);
    void swap(int, int);
    const boost::multi_array<int,2> getIDs();
    const boost::multi_array<int,2> getOcc();
    
    void printMap();

    void setOccupation(boost::multi_array<int,2>);

};

#endif // CHOLESTERINSYSTEM_H
