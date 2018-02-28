#include "../lib/enhance.hpp"
#include "system/mchost.h"
#include "system/datafile.h"
#include "system/omegaoptimizer.h"
#include "system/inputfile.h"

#include <random>
#include <iostream>

int main(int argc, char **argv) {

    #ifndef NDEBUG
    std::cout<<"Debugging"<<std::endl;
    #endif
    

    enhance::seed = std::random_device{}();
//     enhance::seed = 123456749;
    enhance::rand_engine.seed(enhance::seed);

    InputFile::loadFile("in.txt");
    
    
    MCHost mchost;
    mchost.setup();
    mchost.run();

//     OmegaOptimizer omegaoptimizer;
//     omegaoptimizer.setupOptimization();
//     omegaoptimizer.optimizeOmega();
    


    
    
    
    
    return 0;
}
