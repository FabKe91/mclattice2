#include "mchost.h"

void MCHost::setup()
{
    #ifndef NDEBUG
    std::cout<<"MCHost::setup"<<std::endl;
    #endif
    
    
    lipidproperties.reset(new LipidProperties());
    lipidproperties->readParas();
    

    lipidsystem.readParas(lipidproperties);
    lipidsystem.setup();
    
    
    datafile.reset(new DataFile(lipidsystem,InputFile::outs));
    datafile->createFile();   
    

    



    steps=InputFile::paras.at("steps");
    imageRate=InputFile::paras.at("imageRate");
}



void MCHost::run()
{
    #ifndef NDEBUG
    std::cout<<"MCHost::run"<<std::endl;
    #endif
    
    auto start = std::chrono::system_clock::now();
    unsigned int t=0;
//     int loopCounter=0;


    //     while(t<=steps)

    for(int loopCounter=0;loopCounter<=steps;loopCounter++)
    {   
//         loopCounter++;
        doSystemloop();
        t+=InputFile::paras.at("width")*InputFile::paras.at("height");
        
/*        if(loopCounter==1000)
        {
            std::cout<<"update T!"<<std::endl;
            InputFile::paras["T"]=320;
            InputFile::paras["kBT"]=InputFile::paras.at("T")*InputFile::paras.at("kB");
            lipidproperties->updateKBT();
        } */       
        
        if(loopCounter%imageRate==0)
        {
            datafile->writeStep();
            std::cout<<"loop: "<<loopCounter<<" flucAccepRate: "<< (t-notAcceptedFlucs)/(double)t<<" swapAccepRate: "<< (t-notAcceptedFlucs-notAcceptedSwaps)/(double)(t-notAcceptedFlucs)<<std::endl;
            
            std::chrono::duration<double> elapsed_seconds=std::chrono::system_clock::now()-start;
            
            std::cout<<"time per step: "<<elapsed_seconds.count()/t<<"    mean order: "<<lipidsystem.getMeanOrder()<<std::endl;
        }
    }
    
}

void MCHost::doSystemloop() //loop one time over all lipids
{
    for(unsigned int i=0;i<InputFile::paras.at("width");i++)
    for(unsigned int j=0;j<InputFile::paras.at("height");j++)
    {
        lipidsystem.setHost(i,j);
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
    std::cout<<"MCHost::acceptance:   DeltaG: "<<(FreeEnergie2-FreeEnergie1)<<" DeltaG/kbT: " <<(FreeEnergie2-FreeEnergie1)/InputFile::paras.at("kBT")<<std::endl;
    #endif

    return enhance::random_double(0.0, 1.0) < enhance::fastExp((FreeEnergie1-FreeEnergie2)/InputFile::paras.at("kBT")) ? true : false;
    
}
