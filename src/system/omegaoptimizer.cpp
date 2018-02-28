#include "omegaoptimizer.h"

void OmegaOptimizer::setupOptimization()
{


    
    lipidproperties.reset(new LipidProperties());
    lipidproperties->readParas();
    

    lipidsystem.readParas(lipidproperties);
    lipidsystem.setup();
    
    
    T=InputFile::paras.at("T");
    kB=InputFile::paras.at("kB");
    steps=InputFile::paras.at("steps");
    kBT=kB*T;
    imageRate=InputFile::paras.at("imageRate");
    
    
    
    //change T
//     InputFile::paras["T"]=330;
//     InputFile::paras["kBT"]=InputFile::paras.at("T")*InputFile::paras.at("kB");
//     lipidproperties->updateKBT();
    
    //setting up the md destr
    std::vector<double> coeff={-0.9767356, 8.69286553, -12.7808724, 12.12000201, -21.41776641, 7.14478559};
    
    for(double order=InputFile::paras.at("minOrder");order<InputFile::paras.at("maxOrder")+InputFile::paras.at("DeltaOrder");order+=InputFile::paras.at("DeltaOrder"))
    {
        MDOrderDestr.push_back(std::exp(enhance::polynom(coeff,order)));       
    }
    
    //normalization
    double normSum=0;
    for(auto& Porder: MDOrderDestr)   normSum+=Porder;
    normSum*=InputFile::paras.at("DeltaOrder");
    for(auto& Porder: MDOrderDestr)   Porder/=normSum;
   
}



void OmegaOptimizer::optimizeOmega()
{
    int type=0;

//     initual guess
//     for(int i=0;i<=(int)InputFile::paras["maxOrderIndex"];i++)
//         lipidproperties->entropyFunction[type][i]=std::log(MDOrderDestr[i])+(lipidproperties->enthalpyFunction[type][type][i]*lipidproperties->neighbourFunction[type][i]/2+lipidproperties->selfEnergieFunction[type][i])/kBT;
//     
    //create files for output and close again, not holding files open, to analyse live
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
    int orderCalcRuns=1000;  //number of runs to calc OrderDestr, getting increase 
    double lastMDDiff=INFINITY;
    double alphaFaktor=1; //faktor to decrease alpha 
    while(true)
    {
        runUntilEquilibrium();
        calcCurrentOrderDestr(orderCalcRuns);
        double StepDiff=0; //change per step
        double MDDiff=0; // diff to MD
        OmegaOut.open("OptimzeOut.txt", std::ios_base::app);
        DestrOut.open("DestrOut.txt", std::ios_base::app);

        for(int i=0;i<=(int)InputFile::paras.at("maxOrderIndex");i++) 
        {
            OmegaOut<<lipidproperties->entropyFunction[type][i]<<" "; 
            DestrOut<<currentOrderDestr[i]<<" ";
            
            MDDiff+=(MDOrderDestr[i]-currentOrderDestr[i])*(MDOrderDestr[i]-currentOrderDestr[i]); //diff to md squared

            if (currentOrderDestr[i]!=0) // if ==0 zero division
            { 
                DeltaEnthr=alpha*std::log(MDOrderDestr[i]/currentOrderDestr[i]);
                StepDiff+=DeltaEnthr*DeltaEnthr;
                lipidproperties->entropyFunction[type][i]+=DeltaEnthr;
            }
        }
        
        //following 2 loops only because zero division error from above, making entropyFunction constant in that case
        for(int i=1;i<=(int)InputFile::paras.at("maxOrderIndex");i++) 
        {
            if (currentOrderDestr[i]==0)
            { 
                lipidproperties->entropyFunction[type][i]= lipidproperties->entropyFunction[type][i-1];                           
            }  
        }
        for(int i=(int)InputFile::paras.at("maxOrderIndex")-1;i>=0;i--) 
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
//         std::cout<<"mean order: "<<lipidsystem.getMeanOrder()<<std::endl;
//         if (loopCounter %50==0)std::cout<<lipidsystem.getMeanOrder()<<std::endl;
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
    currentOrderDestr=std::vector<double>(InputFile::paras.at("maxOrderIndex")+1);//resetting and setting to first step
    std::vector<int> thisLoopOrderDestr;
    
    for(int t=0;t<orderCalcRuns;t++) 
    {
        doSystemloop();
        thisLoopOrderDestr=lipidsystem.getOrderDestr();
        for(int i=0;i<=(int)InputFile::paras.at("maxOrderIndex");i++)  currentOrderDestr[i]+=thisLoopOrderDestr[i];
    }
    for(int i=0;i<=(int)InputFile::paras.at("maxOrderIndex");i++)    //normalization
    {
        currentOrderDestr[i]/=InputFile::paras.at("DeltaOrder")*InputFile::paras.at("width")*InputFile::paras.at("height")*orderCalcRuns;       
    }
}

void OmegaOptimizer::doSystemloop() //loop one time over all lipids
{
    for(unsigned int id=0;id<InputFile::paras.at("width")*InputFile::paras.at("height");id++)
    {
        lipidsystem.setHost(id);
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
    std::cout<<"OmegaOptimizer::acceptance:   DeltaG: "<<(FreeEnergie2-FreeEnergie1)<<" DeltaG/kbT: " <<(FreeEnergie2-FreeEnergie1)/kBT<<std::endl;
    #endif

    return enhance::random_double(0.0, 1.0) < enhance::fastExp((FreeEnergie1-FreeEnergie2)/kBT) ? true : false;
    
}



// void OmegaOptimizer::fitOrderDestr() //copy pasta from https://stackoverflow.com/questions/45734608/least-squares-polynomial-fitting-works-only-with-even-number-of-coordinates
// {
// using namespace std;
//     int d, i, j, k, n;
//     d=6; //degree
//     double **A = new double *[d+1];
//     for (k = 0; k < d+1; k++) {
//         A[k] = new double[d+2];
//         for (i = 0; i < d+2; i++) {
//             A[k][i] = 0.0;}}
//     n=(int)InputFile::paras["maxOrderIndex"]+1;
//     double * x = new double[n]; 
//     double * y = new double[n];
//     for (i = 0; i < n; i++)
//         x[i]=InputFile::paras["minOrder"]+InputFile::paras["DeltaOrder"]*i;
//     for (i = 0; i < n; i++)
//         y[i]=logWithZero(currentOrderDestr[i]);
//     for(k = 0; k < d+1; k++){
//         for (i = 0; i < n; i++) {
//             for (j = 0; j < d+1; j++) {
//                 A[k][j] += pow(x[i], j + k);}
//             A[k][d+1] += y[i]*pow(x[i], k);}}
//     for(k = 0; k < d+1; k++){       // invert matrix
//         double q = A[k][k];         //   divide A[k][] by A[k][k]
//                                     //   if q == 0, would need to swap rows                                 
//         for(i = 0; i < d+2; i++){
//             A[k][i] /= q;}
//         for(j = 0; j < d+1; j++){   //   zero out column A[][k]
//             if(j == k)
//                 continue;
//             double m = A[j][k];         
//             for(i = 0; i < d+2; i++){
//                 A[j][i] -= m*A[k][i];}}}
// 
//                 
//     fittetCurrentOrderDestr=std::vector<double>(0);
//     
//     for(double order=InputFile::paras["minOrder"];order<InputFile::paras["maxOrder"]+InputFile::paras["DeltaOrder"];order+=InputFile::paras["DeltaOrder"])
//     {
//         double currentOrder=0;
//         double orderPotenz=1;
//         for(k = 0; k <= d; k++)
//         {
//             currentOrder+=orderPotenz*A[k][d+1];
//             orderPotenz*=order;
//         }
//         fittetCurrentOrderDestr.push_back(std::exp(currentOrder));
//     }
//     
//     for (k = 0; k < d+1; k++)
//         delete[] A[k];
//     delete[]A;
//     delete[]y;
//     delete[]x;
// }

