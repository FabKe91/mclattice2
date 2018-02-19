#include "mchost.h"

void MCHost::setup(std::string filename)
{

    inputfile.reset(new InputFile(filename));

    
    lipidproperties.reset(new LipidProperties());
    lipidproperties->readParas(inputfile);
    

    lipidsystem.readParas(inputfile,lipidproperties);
    lipidsystem.setup();
    
    
    datafile.reset(new DataFile(lipidsystem,inputfile->outs));
    datafile->createFile();   
    

    


    T=inputfile->paras["T"];
    kB=inputfile->paras["kB"];
    steps=inputfile->paras["steps"];
    kBT=kB*T;
    imageRate=inputfile->paras["imageRate"];
    
   
    

}

void MCHost::optimizeOmega()
{
    auto start = std::chrono::system_clock::now();
    unsigned int t=0;
    
    
    
    int meanOrder=0; //for convergence check
    int loopCounter=0;
    
    while(true)
    {
        for(unsigned int i=0;i<inputfile->paras["width"];i++)
        for(unsigned int j=0;j<inputfile->paras["height"];j++)
        {
            t++;
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
            
            
            else
            {
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
            
            if(t%imageRate==0)
            {
            datafile->writeStep(lipidsystem);
            std::cout<<"t: "<<t<<" flucAccepRate: "<< (t-notAcceptedFlucs)/(double)t<<" swapAccepRate: "<< (t-notAcceptedFlucs-notAcceptedSwaps)/(double)(t-notAcceptedFlucs)<<std::endl;
            std::chrono::duration<double> elapsed_seconds=std::chrono::system_clock::now()-start;
            std::cout<<"time per step: "<<elapsed_seconds.count()/t<<"    mean order: "<<lipidsystem.getMeanOrder(0)<<std::endl;
            }
        }
    }
}


void MCHost::run()
{
    auto start = std::chrono::system_clock::now();
    unsigned int t=0;
    
    while(t<=steps)
    for(unsigned int i=0;i<inputfile->paras["width"] and t<=steps;i++)
    for(unsigned int j=0;j<inputfile->paras["height"] and t<=steps;j++)
    {
        t++;
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
        
        
        else
        {
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
        
        if(t%imageRate==0)
        {
        datafile->writeStep(lipidsystem);
        std::cout<<"t: "<<t<<" flucAccepRate: "<< (t-notAcceptedFlucs)/(double)t<<" swapAccepRate: "<< (t-notAcceptedFlucs-notAcceptedSwaps)/(double)(t-notAcceptedFlucs)<<std::endl;
        std::chrono::duration<double> elapsed_seconds=std::chrono::system_clock::now()-start;
        std::cout<<"time per step: "<<elapsed_seconds.count()/t<<"    mean order: "<<lipidsystem.getMeanOrder(0)<<std::endl;
        }
    }
}





bool MCHost::acceptance(const double FreeEnergie1, const double FreeEnergie2)
{
    #ifndef NDEBUG
    std::cout<<"MCHost::acceptance:   DeltaG: "<<(FreeEnergie2-FreeEnergie1)<<" DeltaG/kbT: " <<(FreeEnergie2-FreeEnergie1)/kBT<<std::endl;
    #endif

    return enhance::random_double(0.0, 1.0) < enhance::fastExp((FreeEnergie1-FreeEnergie2)/kBT) ? true : false;
    
}
