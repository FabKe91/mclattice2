#include "omegaoptimizer.h"


    #define ID0 lastSwappedIDs[0]
    #define LIPID0 lipidsystem.lipids[ID0]
    #define posX0 LIPID0.posX
    #define posY0 LIPID0.posY
    
    #define ID1 lastSwappedIDs[1]
    #define LIPID1 lipidsystem.lipids[ID1]
    #define posX1 LIPID1.posX
    #define posY1 LIPID1.posY
    

void OmegaOptimizer::setupOptimization(std::string inputFileName,std::string typeName)
{
    #ifndef NDEBUG
    std::cout<<"OmegaOptimizer::setupOptimization"<<std::endl;
    #endif
    
    inputfile.reset(new InputFile(inputFileName));
    
    
    lipidproperties.reset(new LipidProperties());
    lipidproperties->readParas(inputfile);
    
    type=inputfile->typeMap[typeName]; //get type number from typeName
    for (auto& conc: inputfile->concentrations) conc=0;    //setting all conc to 0
    inputfile->concentrations[type]=1;     //setting conc of type to 1


    

    lipidsystem.readParas(lipidproperties,inputfile);
    lipidsystem.setup();
    
    
    for(int i=0;i<inputfile->width*inputfile->height;i++) IDs.push_back(i); //setup IDs to get rndOrder in each doSystemloop
    
    

    //setting up the md distr
    std::vector<double> coeff;
    if (typeName=="DPPC")   coeff=std::vector<double>({-0.9767356, 8.69286553, -12.7808724, 12.12000201, -21.41776641, 7.14478559});
    else if (typeName=="DUPC") coeff=std::vector<double>({-0.14122, 7.51277, -9.36903, -4.43679, -97.86418, 192.92704, 19.37517, -168.20577});
    else throw std::invalid_argument( "no MD distr found for type: "+typeName);
        
    for(double order=inputfile->paras.at("minOrder");order<inputfile->paras.at("maxOrder")+inputfile->paras.at("DeltaOrder");order+=inputfile->paras.at("DeltaOrder"))
    {
        MDOrderDistr.push_back(std::exp(enhance::polynom(coeff,order)));       
    }
    
    //normalization
    double normSum=0;
    for(int i=inputfile->types[type].minOrder;i<=inputfile->types[type].maxOrder;i++)   normSum+=MDOrderDistr[i];
    normSum*=inputfile->paras.at("DeltaOrder");
    for(auto& Porder: MDOrderDistr)   Porder/=normSum;
   
}



void OmegaOptimizer::optimizeOmega()
{

//     initual guess
    for(int i=0;i<=(int)inputfile->paras["maxOrderIndex"];i++)
        lipidproperties->entropyFunction[type][i]=std::log(MDOrderDistr[i])+(lipidproperties->enthalpyFunction[type][type][i]*lipidproperties->neighbourFunction[type][i]/2+lipidproperties->selfEnergieFunction[type][i])/inputfile->kBT;
    
    //create files for output
    std::ofstream OmegaOut;
    OmegaOut.open("OptimzeOut.txt", std::ios_base::out);
    std::ofstream DistrOut;
    DistrOut.open("DistrOut.txt", std::ios_base::out);
    OmegaOut.close();
    DistrOut.close();
    
    int run=0;
    double DeltaEnthr=0; 
    double maxAlpha=0.5;
    double alpha=maxAlpha;
    int orderCalcRuns=1000;  //number of runs to calc OrderDistr
    double lastMDDiff=INFINITY;
    
    while(true)
    {
        runUntilEquilibrium();
        calcCurrentOrderDistr(orderCalcRuns);
        double StepDiff=0; //diff to last step
        double MDDiff=0; // diff to MD
        OmegaOut.open("OptimzeOut.txt", std::ios_base::app);
        DistrOut.open("DistrOut.txt", std::ios_base::app);

        int VZW=0;
        bool VZ=true; //true=positiv
        bool lastVZ=true;
        for(int i=0;i<=(int)inputfile->paras.at("maxOrderIndex");i++) 
        {
            //write to file
            OmegaOut<<lipidproperties->entropyFunction[type][i]<<" "; 
            DistrOut<<currentOrderDistr[i]<<" ";
            
            MDDiff+=(MDOrderDistr[i]-currentOrderDistr[i])*(MDOrderDistr[i]-currentOrderDistr[i]); //diff to md squared

            //update omega
            if (currentOrderDistr[i]!=0) // if ==0 zero division
            { 
                DeltaEnthr=alpha*std::log(MDOrderDistr[i]/currentOrderDistr[i]);
                VZ=DeltaEnthr>0;
                if(VZ != lastVZ)
                {
                    VZW++;
                    lastVZ=VZ;
                }
                StepDiff+=DeltaEnthr*DeltaEnthr;
                lipidproperties->entropyFunction[type][i]+=DeltaEnthr;
            }
        }
        
        //following 2 loops only because zero division error from above, making entropyFunction constant in that case
        for(int i=1;i<=(int)inputfile->paras.at("maxOrderIndex");i++) 
        {
            if (currentOrderDistr[i]==0)
            { 
                lipidproperties->entropyFunction[type][i]= lipidproperties->entropyFunction[type][i-1];                           
            }
        }
        for(int i=(int)inputfile->paras.at("maxOrderIndex")-1;i>=0;i--) 
        {
            if (currentOrderDistr[i]==0)
            { 
                lipidproperties->entropyFunction[type][i]= lipidproperties->entropyFunction[type][i+1];                           
            }
        }

        OmegaOut<<"\n";
        DistrOut<<"\n";
        OmegaOut.close();
        DistrOut.close();

        if (MDDiff>lastMDDiff)
            {
                orderCalcRuns*=2;
                if(VZW <15)
                    alpha*=0.5;
            }
        if (alpha > maxAlpha) alpha=maxAlpha;
        std::cout<<"run "<<run<<" StepDiff "<<StepDiff<<" MDDiff "<<MDDiff<<" orderCalcRuns "<<orderCalcRuns<<" alpha "<<alpha<<" VZW "<<VZW<<std::endl;

        lastMDDiff=MDDiff;
        run++;

    }
}

void OmegaOptimizer::runUntilEquilibrium()
{ 
    int meanOrder=0; //for convergence check
    int loopCounter=0;
    
    while(true) //run until equi
    {
        doSystemloop();
        loopCounter++;
        if (loopCounter %100==0)
        {
            std::cout<<"loops: "<<loopCounter<<" mean order:  last: "<<meanOrder<<" now: "<<lipidsystem.getMeanOrder()<<std::endl;
            
            if (std::abs(lipidsystem.getMeanOrder()-meanOrder)<=1)//check convergence
            {
                std::cout<<"Equilibrium!"<<std::endl;
                break;
            }
            meanOrder=lipidsystem.getMeanOrder();

        }
    }
}

void OmegaOptimizer::calcCurrentOrderDistr(int orderCalcRuns)//doing orderCalcRuns more loops to calc order distr
{
    currentOrderDistr=std::vector<double>(inputfile->paras.at("maxOrderIndex")+1);//resetting and setting to first step
    std::vector<int> thisLoopOrderDistr;
    
    for(int t=0;t<orderCalcRuns;t++) 
    {
        doSystemloop();
        thisLoopOrderDistr=lipidsystem.getOrderDistr();
        for(int i=0;i<=(int)inputfile->paras.at("maxOrderIndex");i++)  currentOrderDistr[i]+=thisLoopOrderDistr[i];
    }
    for(int i=0;i<=(int)inputfile->paras.at("maxOrderIndex");i++)    //normalization
    {
        currentOrderDistr[i]/=inputfile->paras.at("DeltaOrder")*inputfile->paras.at("width")*inputfile->paras.at("height")*orderCalcRuns;       
    }
}

void OmegaOptimizer::doSystemloop() //loop one time over all lipids
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

double OmegaOptimizer::calcSwapEnthalpy()
{
    #ifndef NDEBUG
    std::cout<<"OmegaOptimizer::calcSwapEnthalpy"<<std::endl;
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

double OmegaOptimizer::calcHostFreeEnerg()
{
    #ifndef NDEBUG
    std::cout<<"OmegaOptimizer::calcHostFreeEnerg"<<std::endl;
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


void OmegaOptimizer::setHost(int x, int y)
{
    ID0=lipidsystem.map[x][y];
}

void OmegaOptimizer::setHost(int ID)
{
    ID0=ID;
}

void OmegaOptimizer::setRNDHost()
{
    ID0=enhance::random_int(0,inputfile->height*inputfile->width-1);
}

void OmegaOptimizer::setPartner()
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


bool OmegaOptimizer::acceptance(const double FreeEnergie1, const double FreeEnergie2)
{
    #ifndef NDEBUG
    std::cout<<"OmegaOptimizer::acceptance:   DeltaG: "<<(FreeEnergie2-FreeEnergie1)<<" DeltaG/kbT: " <<(FreeEnergie2-FreeEnergie1)/inputfile->kBT<<std::endl;
    #endif

    return enhance::random_double(0.0, 1.0) < enhance::fastExp((FreeEnergie1-FreeEnergie2)/inputfile->kBT);
    
}

