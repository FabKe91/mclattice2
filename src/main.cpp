#include "system/mchost.h"
#include "system/datafile.h"

#include "../lib/enhance.hpp"
#include <random>
#include <iostream>

int main(int argc, char **argv) {

    #ifndef NDEBUG
    std::cout<<"Debugging"<<std::endl;
    #endif
    

    MCHost mchost;
    //  enhance::seed = std::random_device{}();
    enhance::seed = 123456789;
    enhance::rand_engine.seed(enhance::seed);

    mchost.setup("in.txt");
//     mchost.run();
    mchost.optimizeOmega();
    


    
    
    
    
    return 0;
}
