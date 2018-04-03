#include "omegaoptimizer.h"


    #define HOST_LIPID lipidsystem.lipids[hostID]
    #define HOST_posX HOST_LIPID.posX
    #define HOST_posY HOST_LIPID.posY
    
    #define PARTNER_LIPID lipidsystem.lipids[partnerID]
    #define PARTNER_posX PARTNER_LIPID.posX
    #define PARTNER_posY PARTNER_LIPID.posY
    

void OmegaOptimizer::setupOptimization(std::string typeName)
{
    #ifndef NDEBUG
    std::cout<<"OmegaOptimizer::setupOptimization"<<std::endl;
    #endif
    
    inputfile.reset(new InputFile());
    inputfile->T=330;
    inputfile->paras["T"]=330;
    inputfile->kBT=inputfile->T*inputfile->paras.at("kB");
    
    
    lipidproperties.reset(new LipidProperties());
    lipidproperties->readParas(inputfile);
    
    type=inputfile->typeMap[typeName]; //get type number from typeName
    for (auto& conc: inputfile->concentrations) conc=0;    //setting all conc to 0
    inputfile->concentrations[type]=1;     //setting conc of type to 1


    

    lipidsystem.readParas(lipidproperties,inputfile);
    lipidsystem.setup();
    
    cholesterinsystem.setup(inputfile);

    
    for(int i=0;i<inputfile->width*inputfile->height;i++)
    {
        IDs.push_back(i);
        cholIDs.push_back(i);
    }
    
    

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
    for(int i=std::get<1>(inputfile->types[type]);i<=std::get<2>(inputfile->types[type]);i++)   normSum+=MDOrderDistr[i];
    normSum*=inputfile->paras.at("DeltaOrder");
    for(auto& Porder: MDOrderDistr)   Porder/=normSum;
   
}



void OmegaOptimizer::optimizeOmega()
{
    #ifndef NDEBUG
    std::cout<<"OmegaOptimizer::optimizeOmega"<<std::endl;
    #endif
    
//     initual guess
    for(int i=0;i<=(int)inputfile->paras["maxOrderIndex"];i++)
        lipidproperties->entropyFunction[type][i]=std::log(MDOrderDistr[i])+(lipidproperties->enthalpyFunction[type][type][0][i]*lipidproperties->neighbourFunction[type][i]/2+lipidproperties->selfEnergieFunction[type][i])/inputfile->kBT; //die [type][type][0] ist eig nicht sinnvoll, da nicht ueberall 0 nachbarn
    
    //create files for output
    std::ofstream OmegaOut;
    OmegaOut.open("OptimzeOut.txt", std::ios_base::out);
    std::ofstream DistrOut;
    DistrOut.open("DistrOut.txt", std::ios_base::out);
    OmegaOut.close();
    DistrOut.close();
    
    int run=0;
    double DeltaEnthr=0; 
    double alpha=0.2;
    int orderCalcRuns=200;  //number of runs to calc OrderDistr
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
            }
//         if (alpha > maxAlpha) alpha=maxAlpha;
        alpha*=std::pow(MDDiff/StepDiff/2,0.5);
        std::cout<<"run "<<run<<" StepDiff "<<StepDiff<<" MDDiff "<<MDDiff<<" orderCalcRuns "<<orderCalcRuns<<" alpha "<<alpha<<" VZW "<<VZW<<std::endl;

        lastMDDiff=MDDiff;
        run++;

    }
}

void OmegaOptimizer::runUntilEquilibrium()
{ 
    #ifndef NDEBUG
    std::cout<<"OmegaOptimizer::runUntilEquilibrium"<<std::endl;
    #endif
    
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

