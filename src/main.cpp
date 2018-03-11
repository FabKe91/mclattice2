#include "../lib/enhance.hpp"
#include "system/mchost.h"
#include "system/datafile.h"
#include "system/omegaoptimizer.h"
#include "system/inputfile.h"

#include <random>
#include <iostream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char **argv) {

    
    
    enhance::seed = std::random_device{}();
//     enhance::seed = 123456749;

    #ifndef NDEBUG
    std::cout<<"Debugging"<<std::endl;
    enhance::seed = 123456749;
    #endif
    
    enhance::rand_engine.seed(enhance::seed);

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("optimize", po::value<std::string>(), "do optimization, type name")
        ("o", po::value<std::string>(), "do optimization, type name")
        ("continue","continue from existing file named 'out.h5'")
    ;
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    
    
    
    
    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        return 1;
    }
    else if (vm.count("optimize") or vm.count("o"))
    {
        std::string typeName;
        if (vm.count("optimize")) typeName=vm["optimize"].as<std::string>();
        if (vm.count("o")) typeName=vm["o"].as<std::string>();

        OmegaOptimizer omegaoptimizer;
        omegaoptimizer.setupOptimization("in.txt",typeName);
        omegaoptimizer.optimizeOmega();
    }
    else if (vm.count("continue"))
    {
        MCHost mchost;
        mchost.setupForRestart("in.txt");
        mchost.run();
    }
    else
    {
        MCHost mchost;
        mchost.setup("in.txt");
        mchost.run();
    }
  
    
    
    
    return 0;
}
