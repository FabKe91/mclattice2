#include "mchost.h"

void MCHost::setup(std::string inputFileName)
{
    #ifndef NDEBUG
    std::cerr<<"MCHost::setup"<<std::endl;
    #endif
    
    
    inputfile.reset(new InputFile(inputFileName));  //all input parameters are stored in the shared pointer "inputfile". all classes get the pointer 
    
    lipidproperties.reset(new LipidProperties()); //the fit functions are stored in lipidproperties
    lipidproperties->readParas(inputfile);
    

    lipidsystem.readParas(lipidproperties,inputfile); 
    lipidsystem.setup();
    
    
    datafile.reset(new DataFile(lipidsystem,inputfile));
    datafile->createFile();   
    

    


    for(int i=0;i<inputfile->width*inputfile->height;i++) IDs.push_back(i);
    
    steps=inputfile->paras.at("steps");
    imageRate=inputfile->paras.at("imageRate");
}



void MCHost::run()
{
    #ifndef NDEBUG
    std::cerr<<"MCHost::run"<<std::endl;
    #endif
    
    auto startTime = std::chrono::system_clock::now();
    auto lastTime = std::chrono::system_clock::now();
    long long t=0;
    
    
    
    datafile->writeStep();
    for(int loopCounter=1;loopCounter<=steps;loopCounter++)
    {   
        doSystemloop();
        t+=inputfile->paras.at("width")*inputfile->paras.at("height");
        
/*        if(loopCounter==1000)
        {
            std::cout<<"update T!"<<std::endl;
            inputfile->paras["T"]=320;
            inputfile->paras["kBT"]=inputfile->paras.at("T")*inputfile->paras.at("kB");
            lipidproperties->updateKBT();
        } */       
        
        if(loopCounter%imageRate==0)
        {
            datafile->writeStep();
            std::cout<<"loop: "<<loopCounter<<" flucAccepRate: "<< (t-notAcceptedFlucs)/(double)t<<" swapAccepRate: "<< (t-notAcceptedFlucs-notAcceptedSwaps)/(double)(t-notAcceptedFlucs)<<std::endl;
            
            auto currTime = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds_start=currTime-startTime;
            std::chrono::duration<double> elapsed_seconds_last=currTime-lastTime;
            std::time_t end_time = std::chrono::system_clock::to_time_t(currTime);

            
            std::cout<<"time per step curr: "<<elapsed_seconds_last.count()/imageRate/inputfile->paras.at("width")/inputfile->paras.at("height")<<" mean: "<<elapsed_seconds_start.count()/t<<"    mean order: "<<lipidsystem.getMeanOrder()<<" "<<std::ctime(&end_time)<<std::endl;
            lastTime = currTime;

        }
    }
    
}

void MCHost::doSystemloop() //loop one time over all lipids
{
    std::shuffle(IDs.begin(), IDs.end(), enhance::rand_engine);
    
    for(int &ID: IDs)
    {
        lipidsystem.setHost(ID);
        double FreeEnergie1=0;
        double FreeEnergie2=0;

        FreeEnergie1=lipidsystem.calcHostFreeEnerg();
        lipidsystem.fluctuate();
        FreeEnergie2=lipidsystem.calcHostFreeEnerg();

        if(!acceptance(FreeEnergie1,FreeEnergie2))
        {
            lipidsystem.fluctuateBack();
            notAcceptedFlucs++;

        }

        lipidsystem.setPartner();
        FreeEnergie1=lipidsystem.calcSwapEnthalpy();
        lipidsystem.swap();
        FreeEnergie2=lipidsystem.calcSwapEnthalpy();

        if(!acceptance(FreeEnergie1,FreeEnergie2))
        {
            lipidsystem.swap();
            notAcceptedSwaps++;

        }
    }
}




bool MCHost::acceptance(const double FreeEnergie1, const double FreeEnergie2)
{
    #ifndef NDEBUG
    std::cerr<<"MCHost::acceptance:   DeltaG: "<<(FreeEnergie2-FreeEnergie1)<<" DeltaG/kbT: " <<(FreeEnergie2-FreeEnergie1)/inputfile->kBT<<std::endl;
    #endif

    return enhance::random_double(0.0, 1.0) < enhance::fastExp((FreeEnergie1-FreeEnergie2)/inputfile->kBT) ? true : false;
    
}
