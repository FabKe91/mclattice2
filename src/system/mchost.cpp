#include "mchost.h"


    #define ID0 lastSwappedIDs[0]
    #define LIPID0 lipidsystem.lipids[ID0]
    #define posX0 LIPID0.posX
    #define posY0 LIPID0.posY
    
    #define ID1 lastSwappedIDs[1]
    #define LIPID1 lipidsystem.lipids[ID1]
    #define posX1 LIPID1.posX
    #define posY1 LIPID1.posY
    

    
void MCHost::setup(std::string inputFileName)
{
    #ifndef NDEBUG
    std::cout<<"MCHost::setup"<<std::endl;
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

void MCHost::setupForRestart(std::string inputFileName)
{
    #ifndef NDEBUG
    std::cout<<"MCHost::setup"<<std::endl;
    #endif
    
    
    inputfile.reset(new InputFile(inputFileName));  //all input parameters are stored in the shared pointer "inputfile". all classes get the pointer 
    
    
    datafile.reset(new DataFile(lipidsystem,inputfile));
    datafile->readFile();
    
    if(inputfile->paras.at("steps")/inputfile->paras.at("imageRate")+1==datafile->getImages()) throw std::logic_error("run already finished");
    loopCounter=(datafile->getImages()-1)*inputfile->paras.at("imageRate");
    

    lipidproperties.reset(new LipidProperties()); //the fit functions are stored in lipidproperties
    lipidproperties->readParas(inputfile);
    

    lipidsystem.readParas(lipidproperties,inputfile); 
    lipidsystem.setup();
    
    
    lipidsystem.setOrder(datafile->getLastStep("orderPara"));
    lipidsystem.setTypes(datafile->getLastStep("Type"));
//     lipidsystem.setIDs(datafile->getLastStep("IDs"));



    for(int i=0;i<inputfile->width*inputfile->height;i++) IDs.push_back(i);
    
    steps=inputfile->paras.at("steps");
    imageRate=inputfile->paras.at("imageRate");
}


void MCHost::run()
{
    #ifndef NDEBUG
    std::cout<<"MCHost::run"<<std::endl;
    #endif
    
    auto startTime = std::chrono::system_clock::now();
    auto lastTime = std::chrono::system_clock::now();
    long long t=0;
    
    
    if (loopCounter == 0) datafile->writeStep();
    while(loopCounter<steps)
    {   
        loopCounter++;
        doSystemloop();
        t+=inputfile->width*inputfile->height;

        if(loopCounter%imageRate==0)
        {
            datafile->writeStep();
            
            std::cout<<"loop: "<<loopCounter<<" flucAccepRate: "<< (t-notAcceptedFlucs)/(double)t<<" swapAccepRate: "<< (t-notAcceptedFlucs-notAcceptedSwaps)/(double)(t-notAcceptedFlucs)<<std::endl;
            
            auto currTime = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds_start=currTime-startTime;
            std::chrono::duration<double> elapsed_seconds_last=currTime-lastTime;
            std::time_t end_time = std::chrono::system_clock::to_time_t(currTime);
            
            std::cout<<"time per step curr: "<<elapsed_seconds_last.count()/imageRate/inputfile->width/inputfile->height<<" mean: "<<elapsed_seconds_start.count()/t<<"    mean order: "<<lipidsystem.getMeanOrder()<<" "<<std::ctime(&end_time)<<std::endl;
            lastTime = currTime;
        }
    }
}

void MCHost::doSystemloop() //loop one time over all lipids
{
    std::shuffle(IDs.begin(), IDs.end(), enhance::rand_engine);
    
    for(int &ID: IDs)
    {
        setHost(ID);
        double FreeEnergie1=0;
        double FreeEnergie2=0;

        FreeEnergie1=calcHostFreeEnerg();
        lipidsystem.fluctuate(ID0);
        FreeEnergie2=calcHostFreeEnerg();

        if(!acceptance(FreeEnergie1,FreeEnergie2))
        {
            lipidsystem.fluctuateBack(ID0);
            notAcceptedFlucs++;
        }

        setPartner();
        FreeEnergie1=calcSwapEnthalpy();
        lipidsystem.swap(ID0,ID1);
        FreeEnergie2=calcSwapEnthalpy();

        if(!acceptance(FreeEnergie1,FreeEnergie2))
        {
            lipidsystem.swap(ID0,ID1);
            notAcceptedSwaps++;
        }
    }
}

double MCHost::calcSwapEnthalpy()
{
    #ifndef NDEBUG
    std::cout<<"MCHost::calcSwapEnthalpy"<<std::endl;
    #endif

    double H=0;
    
    switch(rdnPartnerNumber)
    {
        case 0: H+=lipidsystem.calcPairEnthalpy(ID0,lipidsystem.map[(posX0-1+inputfile->width)%inputfile->width][posY0]);
                H+=lipidsystem.calcPairEnthalpy(ID0,lipidsystem.map[posX0][(posY0+1)%inputfile->height]);
                H+=lipidsystem.calcPairEnthalpy(ID0,lipidsystem.map[posX0][(posY0-1+inputfile->height)%inputfile->height]);
                
                H+=lipidsystem.calcPairEnthalpy(ID1,lipidsystem.map[(posX1+1)%inputfile->width][posY1]);
                H+=lipidsystem.calcPairEnthalpy(ID1,lipidsystem.map[posX1][(posY1+1)%inputfile->height]);
                H+=lipidsystem.calcPairEnthalpy(ID1,lipidsystem.map[posX1][(posY1-1+inputfile->height)%inputfile->height]);
                break;
                
        case 1: H+=lipidsystem.calcPairEnthalpy(ID0,lipidsystem.map[(posX0+1)%inputfile->width][posY0]);
                H+=lipidsystem.calcPairEnthalpy(ID0,lipidsystem.map[posX0][(posY0+1)%inputfile->height]);
                H+=lipidsystem.calcPairEnthalpy(ID0,lipidsystem.map[posX0][(posY0-1+inputfile->height)%inputfile->height]);
                
                H+=lipidsystem.calcPairEnthalpy(ID1,lipidsystem.map[(posX1-1+inputfile->width)%inputfile->width][posY1]);
                H+=lipidsystem.calcPairEnthalpy(ID1,lipidsystem.map[posX1][(posY1+1)%inputfile->height]);
                H+=lipidsystem.calcPairEnthalpy(ID1,lipidsystem.map[posX1][(posY1-1+inputfile->height)%inputfile->height]);
                break;
       
        case 2: H+=lipidsystem.calcPairEnthalpy(ID0,lipidsystem.map[(posX0+1)%inputfile->width][posY0]);
                H+=lipidsystem.calcPairEnthalpy(ID0,lipidsystem.map[(posX0-1+inputfile->width)%inputfile->width][posY0]);
                H+=lipidsystem.calcPairEnthalpy(ID0,lipidsystem.map[posX0][(posY0-1+inputfile->height)%inputfile->height]);
                
                H+=lipidsystem.calcPairEnthalpy(ID1,lipidsystem.map[(posX1+1)%inputfile->width][posY1]);
                H+=lipidsystem.calcPairEnthalpy(ID1,lipidsystem.map[(posX1-1+inputfile->width)%inputfile->width][posY1]);
                H+=lipidsystem.calcPairEnthalpy(ID1,lipidsystem.map[posX1][(posY1+1)%inputfile->height]);
                break;
                
        case 3: H+=lipidsystem.calcPairEnthalpy(ID0,lipidsystem.map[(posX0+1)%inputfile->width][posY0]);
                H+=lipidsystem.calcPairEnthalpy(ID0,lipidsystem.map[(posX0-1+inputfile->width)%inputfile->width][posY0]);
                H+=lipidsystem.calcPairEnthalpy(ID0,lipidsystem.map[posX0][(posY0+1)%inputfile->height]);
                
                H+=lipidsystem.calcPairEnthalpy(ID1,lipidsystem.map[(posX1+1)%inputfile->width][posY1]);
                H+=lipidsystem.calcPairEnthalpy(ID1,lipidsystem.map[(posX1-1+inputfile->width)%inputfile->width][posY1]);
                H+=lipidsystem.calcPairEnthalpy(ID1,lipidsystem.map[posX1][(posY1-1+inputfile->height)%inputfile->height]);
    }
    #ifndef NDEBUG
    std::cout<<"H "<<H<<std::endl;
    #endif
    
    return H;
}

double MCHost::calcHostFreeEnerg()
{
    #ifndef NDEBUG
    std::cout<<"MCHost::calcHostFreeEnerg"<<std::endl;
    std::cout<<"ID0 "<<ID0<<" posX0 "<<posX0<<" posY0 "<<posY0<<std::endl;
    #endif

    
    double G=0;

    G+=lipidsystem.calcPairEnthalpy(ID0,lipidsystem.map[(posX0-1+inputfile->width)%inputfile->width][posY0]);
    G+=lipidsystem.calcPairEnthalpy(ID0,lipidsystem.map[(posX0+1)%inputfile->width][posY0]);
    G+=lipidsystem.calcPairEnthalpy(ID0,lipidsystem.map[posX0][(posY0+1)%inputfile->height]);
    G+=lipidsystem.calcPairEnthalpy(ID0,lipidsystem.map[posX0][(posY0-1+inputfile->height)%inputfile->height]);

   
    #ifndef NDEBUG
    std::cout<<"H "<<G<<std::endl;
    std::cout<<"kB T S "<<-inputfile->kBT*lipidproperties->entropyFunction[LIPID0.getType()][LIPID0.getOrder()]<<std::endl;
    std::cout<<"self E "<<lipidproperties->selfEnergieFunction[LIPID0.getType()][LIPID0.getOrder()]<<std::endl;
    #endif

    G+=lipidproperties->selfEnergieFunction[LIPID0.getType()][LIPID0.getOrder()]-inputfile->kBT*lipidproperties->entropyFunction[LIPID0.getType()][LIPID0.getOrder()];
    
    
    return G;
}


void MCHost::setHost(int x, int y)
{
    ID0=lipidsystem.map[x][y];
}

void MCHost::setHost(int ID)
{
    ID0=ID;
}

void MCHost::setRNDHost()
{
    ID0=enhance::random_int(0,inputfile->height*inputfile->width-1);
}

void MCHost::setPartner()
{
    rdnPartnerNumber=enhance::random_int(0,3);
    
    switch (rdnPartnerNumber)
    {
        case 0 :    ID1=lipidsystem.map[(posX0+1)%inputfile->width][posY0]; //up
                    break;
        case 1 :    ID1=lipidsystem.map[(posX0-1+inputfile->width)%inputfile->width][posY0]; //down
                    break;
        case 2 :    ID1=lipidsystem.map[posX0][(posY0+1)%inputfile->height]; //right
                    break;
        case 3 :    ID1=lipidsystem.map[posX0][(posY0-1+inputfile->height)%inputfile->height]; //left
    }
}


bool MCHost::acceptance(const double FreeEnergie1, const double FreeEnergie2)
{
    #ifndef NDEBUG
    std::cout<<"MCHost::acceptance:   DeltaG: "<<(FreeEnergie2-FreeEnergie1)<<" DeltaG/kbT: " <<(FreeEnergie2-FreeEnergie1)/inputfile->kBT<<std::endl;
    #endif

    return enhance::random_double(0.0, 1.0) < enhance::fastExp((FreeEnergie1-FreeEnergie2)/inputfile->kBT) ? true : false;
    
}
