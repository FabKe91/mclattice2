#include "omegaoptimizer.h"

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
    
    

    //setting up the md destr
    std::vector<double> coeff;
    if (typeName=="DPPC")   coeff=std::vector<double>({-0.9767356, 8.69286553, -12.7808724, 12.12000201, -21.41776641, 7.14478559});
    else if (typeName=="DUPC") coeff=std::vector<double>({-0.14122, 7.51277, -9.36903, -4.43679, -97.86418, 192.92704, 19.37517, -168.20577});
    else throw std::invalid_argument( "no MD destr found for type: "+typeName);
        
    for(double order=inputfile->paras.at("minOrder");order<inputfile->paras.at("maxOrder")+inputfile->paras.at("DeltaOrder");order+=inputfile->paras.at("DeltaOrder"))
    {
        MDOrderDestr.push_back(std::exp(enhance::polynom(coeff,order)));       
    }
    
    //normalization
    double normSum=0;
    for(auto& Porder: MDOrderDestr)   normSum+=Porder;
    normSum*=inputfile->paras.at("DeltaOrder");
    for(auto& Porder: MDOrderDestr)   Porder/=normSum;
   
}



void OmegaOptimizer::optimizeOmega()
{

//     initual guess
    for(int i=0;i<=(int)inputfile->paras["maxOrderIndex"];i++)
        lipidproperties->entropyFunction[type][i]=std::log(MDOrderDestr[i])+(lipidproperties->enthalpyFunction[type][type][i]*lipidproperties->neighbourFunction[type][i]/2+lipidproperties->selfEnergieFunction[type][i])/inputfile->kBT;
    
    //create files for output
    std::ofstream OmegaOut;
    OmegaOut.open("OptimzeOut.txt", std::ios_base::out);
    std::ofstream DestrOut;
    DestrOut.open("DestrOut.txt", std::ios_base::out);
    OmegaOut.close();
    DestrOut.close();
    
    int run=0;
    double DeltaEnthr=0; 
    double max_alpha=0.2; 
    double alpha=max_alpha;
    int orderCalcRuns=1000;  //number of runs to calc OrderDestr
    double lastMDDiff=INFINITY;
    double alphaFaktor=1; //faktor to decrease alpha 
    while(true)
    {
        runUntilEquilibrium();
        calcCurrentOrderDestr(orderCalcRuns);
        double StepDiff=0; //diff to last step
        double MDDiff=0; // diff to MD
        OmegaOut.open("OptimzeOut.txt", std::ios_base::app);
        DestrOut.open("DestrOut.txt", std::ios_base::app);

        for(int i=0;i<=(int)inputfile->paras.at("maxOrderIndex");i++) 
        {
            //write to file
            OmegaOut<<lipidproperties->entropyFunction[type][i]<<" "; 
            DestrOut<<currentOrderDestr[i]<<" ";
            
            MDDiff+=(MDOrderDestr[i]-currentOrderDestr[i])*(MDOrderDestr[i]-currentOrderDestr[i]); //diff to md squared

            //update omega
            if (currentOrderDestr[i]!=0) // if ==0 zero division
            { 
                DeltaEnthr=alpha*std::log(MDOrderDestr[i]/currentOrderDestr[i]);
                StepDiff+=DeltaEnthr*DeltaEnthr;
                lipidproperties->entropyFunction[type][i]+=DeltaEnthr;
            }
        }
        
        //following 2 loops only because zero division error from above, making entropyFunction constant in that case
        for(int i=1;i<=(int)inputfile->paras.at("maxOrderIndex");i++) 
        {
            if (currentOrderDestr[i]==0)
            { 
                lipidproperties->entropyFunction[type][i]= lipidproperties->entropyFunction[type][i-1];                           
            }  
        }
        for(int i=(int)inputfile->paras.at("maxOrderIndex")-1;i>=0;i--) 
        {
            if (currentOrderDestr[i]==0)
            { 
                lipidproperties->entropyFunction[type][i]= lipidproperties->entropyFunction[type][i+1];                           
            }  
        }

        OmegaOut<<"\n";
        DestrOut<<"\n";
        OmegaOut.close();
        DestrOut.close();

        if (MDDiff>lastMDDiff) 
        {
            orderCalcRuns*=2;
            alphaFaktor/=2;
            std::cout<<"keine verbesserung! orderCalcRuns = "<<orderCalcRuns<<" alphaFaktor = "<<alphaFaktor<<std::endl;
        }
        alpha=MDDiff*alphaFaktor;
        if (alpha>max_alpha) alpha=max_alpha;
        std::cout<<"run "<<run<<" StepDiff "<<StepDiff<<" MDDiff "<<MDDiff<<" alpha "<<alpha<<std::endl;

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

void OmegaOptimizer::calcCurrentOrderDestr(int orderCalcRuns)//doing orderCalcRuns more loops to calc order destr
{
    currentOrderDestr=std::vector<double>(inputfile->paras.at("maxOrderIndex")+1);//resetting and setting to first step
    std::vector<int> thisLoopOrderDestr;
    
    for(int t=0;t<orderCalcRuns;t++) 
    {
        doSystemloop();
        thisLoopOrderDestr=lipidsystem.getOrderDestr();
        for(int i=0;i<=(int)inputfile->paras.at("maxOrderIndex");i++)  currentOrderDestr[i]+=thisLoopOrderDestr[i];
    }
    for(int i=0;i<=(int)inputfile->paras.at("maxOrderIndex");i++)    //normalization
    {
        currentOrderDestr[i]/=inputfile->paras.at("DeltaOrder")*inputfile->paras.at("width")*inputfile->paras.at("height")*orderCalcRuns;       
    }
}

void OmegaOptimizer::doSystemloop() //loop one time over all lipids
{
    std::shuffle(IDs.begin(), IDs.end(), enhance::rand_engine); //get new rnd id order
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



bool OmegaOptimizer::acceptance(const double FreeEnergie1, const double FreeEnergie2)
{
    #ifndef NDEBUG
    std::cout<<"OmegaOptimizer::acceptance:   DeltaG: "<<(FreeEnergie2-FreeEnergie1)<<" DeltaG/kbT: " <<(FreeEnergie2-FreeEnergie1)/inputfile->kBT<<std::endl;
    #endif

    return enhance::random_double(0.0, 1.0) < enhance::fastExp((FreeEnergie1-FreeEnergie2)/inputfile->kBT);
    
}

