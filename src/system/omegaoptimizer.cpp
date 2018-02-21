#include "omegaoptimizer.h"

void OmegaOptimizer::setupOptimization(std::string filename)
{

    inputfile.reset(new InputFile(filename));

    
    lipidproperties.reset(new LipidProperties());
    lipidproperties->readParas(inputfile);
    

    lipidsystem.readParas(inputfile,lipidproperties);
    lipidsystem.setup();
    
    
    T=inputfile->paras["T"];
    kB=inputfile->paras["kB"];
    steps=inputfile->paras["steps"];
    kBT=kB*T;
    imageRate=inputfile->paras["imageRate"];
    
    
    std::vector<double> coeff={-0.9767356, 8.69286553, -12.7808724, 12.12000201, -21.41776641, 7.14478559};
    for(double order=inputfile->paras["minOrder"];order<inputfile->paras["maxOrder"]+inputfile->paras["DeltaOrder"];order+=inputfile->paras["DeltaOrder"])
    {
        MDOrderDestr.push_back(std::exp(lipidproperties->polynom(coeff,order)));
//         std::cout<<std::exp(lipidproperties->polynom(coeff,order))<<std::endl;
        
    }
}



void OmegaOptimizer::optimizeOmega()
{
    int type=0;

//     initual guess
    for(int i=0;i<=(int)inputfile->paras["maxOrderIndex"];i++)
        lipidproperties->entropyFunction[0][i]=std::log(MDOrderDestr[i])+(lipidproperties->enthalpyFunction[type][type][i]*lipidproperties->neighbourFunction[type][i]/2+lipidproperties->selfEnergieFunction[type][i])/kBT;
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
    int orderCalcRuns=100;
    double lastMDDiff=INFINITY;
    double alphaFaktor=1;
    while(true)
    {
        runUntilEquilibrium();
        calcCurrentOrderDestr(orderCalcRuns);
        double StepDiff=0; //change per step
        double MDDiff=0; // diff to MD
        OmegaOut.open("OptimzeOut.txt", std::ios_base::app);
        DestrOut.open("DestrOut.txt", std::ios_base::app);

        for(int i=0;i<=(int)inputfile->paras["maxOrderIndex"];i++) 
        {
            OmegaOut<<lipidproperties->entropyFunction[type][i]<<" ";
            DestrOut<<currentOrderDestr[i]<<" ";
            MDDiff+=(MDOrderDestr[i]-currentOrderDestr[i])*(MDOrderDestr[i]-currentOrderDestr[i]);

            if (currentOrderDestr[i]!=0)
            { 
                DeltaEnthr=alpha*std::log(MDOrderDestr[i]/currentOrderDestr[i]);
                StepDiff+=DeltaEnthr*DeltaEnthr;
                lipidproperties->entropyFunction[type][i]+=DeltaEnthr;
            }
        }
        for(int i=1;i<=(int)inputfile->paras["maxOrderIndex"];i++) 
        {
            if (currentOrderDestr[i]==0)
            { 
                lipidproperties->entropyFunction[type][i]= lipidproperties->entropyFunction[type][i-1];                           
            }  
        }
        for(int i=(int)inputfile->paras["maxOrderIndex"]-1;i>=0;i--) 
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
                std::cout<<"convergence!"<<std::endl;
                break;
            }
            meanOrder=lipidsystem.getMeanOrder();

        }
    }
}

void OmegaOptimizer::calcCurrentOrderDestr(int orderCalcRuns)//doing orderCalcRuns more loops to calc order destr
{
    currentOrderDestr=std::vector<double>(inputfile->paras["maxOrderIndex"]+1);//resetting and setting to first step
    std::vector<int> thisLoopOrderDestr;
    
    for(int t=0;t<orderCalcRuns;t++) 
    {
        doSystemloop();
        thisLoopOrderDestr=lipidsystem.getOrderDestr();
        for(int i=0;i<=(int)inputfile->paras["maxOrderIndex"];i++)  currentOrderDestr[i]+=thisLoopOrderDestr[i];
    }
    for(int i=0;i<=(int)inputfile->paras["maxOrderIndex"];i++)    //normalization
    {
        currentOrderDestr[i]/=inputfile->paras["DeltaOrder"]*inputfile->paras["width"]*inputfile->paras["height"]*orderCalcRuns;       
    }
}

void OmegaOptimizer::doSystemloop() //loop one time over all lipids
{
    for(unsigned int i=0;i<inputfile->paras["width"];i++)
    for(unsigned int j=0;j<inputfile->paras["height"];j++)
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
            }
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
//     n=(int)inputfile->paras["maxOrderIndex"]+1;
//     double * x = new double[n]; 
//     double * y = new double[n];
//     for (i = 0; i < n; i++)
//         x[i]=inputfile->paras["minOrder"]+inputfile->paras["DeltaOrder"]*i;
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
//     for(double order=inputfile->paras["minOrder"];order<inputfile->paras["maxOrder"]+inputfile->paras["DeltaOrder"];order+=inputfile->paras["DeltaOrder"])
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

