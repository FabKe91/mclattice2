#include "system/mchost.h"
#include "system/datafile.h"
#include "system/omegaoptimizer.h"

#include "../lib/enhance.hpp"
#include <random>
#include <iostream>

int main(int argc, char **argv) {

    #ifndef NDEBUG
    std::cout<<"Debugging"<<std::endl;
    #endif
    

//     enhance::seed = std::random_device{}();
    enhance::seed = 123456789;
    enhance::rand_engine.seed(enhance::seed);

    
    
    MCHost mchost;
    mchost.setup("in.txt");
    mchost.run();

    
//     OmegaOptimizer omegaoptimizer;
//     omegaoptimizer.setupOptimization("in.txt");
//     omegaoptimizer.optimizeOmega();
    


    
    
    
    
    return 0;
}
