#include "../lib/enhance.hpp"
#include "system/mchost.h"
#include "IO/datafile.h"
#include "system/omegaoptimizer.h"
#include "IO/inputfile.h"

#include <random>
#include <iostream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char **argv) {

    
    
    enhance::seed = std::random_device{}();

    #ifndef NDEBUG
    std::cerr<<"Debugging"<<std::endl;
    enhance::seed = 123456749;
    #endif
    
    enhance::rand_engine.seed(enhance::seed);

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("optimize", po::value<std::string>(), "do optimization, type name")
        ("o", po::value<std::string>(), "do optimization, type name")
        ("continue","continue from existing file named 'out.h5'")
        ("normal", "do normal run")
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

        std::cout<<"Starting optimization of "<<typeName<<std::endl;

        OmegaOptimizer omegaoptimizer;
        omegaoptimizer.setupOptimization(typeName);
        omegaoptimizer.optimizeOmega();
    }
    else if (vm.count("continue"))
    {
        MCHost mchost;
        mchost.setupForRestart();
        mchost.run();
    }
    else if ( vm.count("normal") )
    {
        std::cout<<"Starting normal run"<<std::endl;

        MCHost mchost;
        mchost.setup();
        mchost.run();
    }
    else
    {
        std::cout << "Specify run options, or help" << "\n";
    return 1;
    }
  
    
    
    
    return 0;
}
