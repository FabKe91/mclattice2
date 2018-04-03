#include "mchost.h"


    #define HOST_LIPID lipidsystem.lipids[hostID]
    #define HOST_posX HOST_LIPID.posX
    #define HOST_posY HOST_LIPID.posY
    
    #define PARTNER_LIPID lipidsystem.lipids[partnerID]
    #define PARTNER_posX PARTNER_LIPID.posX
    #define PARTNER_posY PARTNER_LIPID.posY
    
    
    #define HOST_CHOL cholesterinsystem.chols[cholHostID]
    #define CHOL_HOST_posX HOST_CHOL.posX
    #define CHOL_HOST_posY HOST_CHOL.posY
    
    #define PARTNER_CHOL cholesterinsystem.chols[cholPartnerID]
    #define CHOL_PARTNER_posX PARTNER_CHOL.posX
    #define CHOL_PARTNER_posY PARTNER_CHOL.posY
    

    
void MCHost::setup()
{
    #ifndef NDEBUG
    std::cout<<"MCHost::setup"<<std::endl;
    #endif
    
    
    inputfile.reset(new InputFile());  //all input parameters are stored in the shared pointer "inputfile". all classes get the pointer 
    
    lipidproperties.reset(new LipidProperties()); //the fit functions are stored in lipidproperties
    lipidproperties->readParas(inputfile);
    

    lipidsystem.readParas(lipidproperties,inputfile); 
    lipidsystem.setup();
    
    cholesterinsystem.setup(inputfile);
    
    
    datafile.reset(new DataFile(lipidsystem, cholesterinsystem, inputfile));
    datafile->createFile();   
    

    


    for(int i=0;i<inputfile->width*inputfile->height;i++)
    {
        IDs.push_back(i);
        cholIDs.push_back(i);
    }
    
    steps=inputfile->paras.at("steps");
    imageRate=inputfile->paras.at("imageRate");
}

void MCHost::setupForRestart()
{
    #ifndef NDEBUG
    std::cout<<"MCHost::setupForRestart"<<std::endl;
    #endif
    
    
    inputfile.reset(new InputFile());  //all input parameters are stored in the shared pointer "inputfile". all classes get the pointer 
    
    
    datafile.reset(new DataFile(lipidsystem, cholesterinsystem, inputfile));
    datafile->readFile();
    
    if(inputfile->paras.at("steps")/inputfile->paras.at("imageRate")+1==datafile->getImages()) throw std::logic_error("run already finished");
    loopCounter=(datafile->getImages()-1)*inputfile->paras.at("imageRate");
    

    lipidproperties.reset(new LipidProperties()); //the fit functions are stored in lipidproperties
    lipidproperties->readParas(inputfile);
    

    lipidsystem.readParas(lipidproperties,inputfile); 
    lipidsystem.setup();
    
    cholesterinsystem.setup(inputfile);

    
    lipidsystem.setOrder(datafile->getLastStep("orderPara"));
    lipidsystem.setTypes(datafile->getLastStep("Type"));
    cholesterinsystem.setOccupation(datafile->getLastStep("Chol"));
//     lipidsystem.setIDs(datafile->getLastStep("IDs"));



    for(int i=0;i<inputfile->width*inputfile->height;i++)
    {
        IDs.push_back(i);
        cholIDs.push_back(i);
    }
    
    
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
            
            std::cout<<"loop: "<<loopCounter<<" flucAccepRate: "<< (t-notAcceptedFlucs)/(double)t<<" swapAccepRate: "<<(t-notAcceptedSwaps)/(double)t<<" chol Hits "<<CholSwaps/(double) t<<" accCholSwaps "<<(CholSwaps-notAcceptedCholSwaps)/(double) CholSwaps<<std::endl;
            
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
    #ifndef NDEBUG
    std::cout<<"MCHost::doSystemloop"<<std::endl;
    #endif
    
    std::shuffle(IDs.begin(), IDs.end(), enhance::rand_engine);
    std::shuffle(cholIDs.begin(), cholIDs.end(), enhance::rand_engine);
    
    for(int i=0;i<inputfile->width*inputfile->height;i++)
    {
        setHost(IDs[i]);
        setPartner();
        double FreeEnergie1=0;
        double FreeEnergie2=0;

        FreeEnergie1=calcHostFreeEnerg();
        lipidsystem.fluctuate(hostID);
        FreeEnergie2=calcHostFreeEnerg();

        if(!acceptance(FreeEnergie1,FreeEnergie2))
        {
            lipidsystem.fluctuateBack(hostID);
            notAcceptedFlucs++;
        }

        FreeEnergie1=calcSwapEnthalpy();
        
        lipidsystem.swap(hostID,partnerID);
        partnerNeighbours=getLipidNeighOfLipid(partnerID);
        hostNeighbours=getLipidNeighOfLipid(hostID);

        
        FreeEnergie2=calcSwapEnthalpy();

        if(!acceptance(FreeEnergie1,FreeEnergie2))
        {
            lipidsystem.swap(hostID,partnerID);
            notAcceptedSwaps++;
        }
        
        
        
        
        if(setCholHost(cholIDs[i]))
        {
            if(!setCholPartner())
            {
                CholSwaps++;
                FreeEnergie1=calcCholSwapEnerg();
                cholesterinsystem.swap(cholHostID,cholPartnerID);
                FreeEnergie2=calcCholSwapEnerg();

                if(!acceptance(FreeEnergie1,FreeEnergie2))
                {
                    cholesterinsystem.swap(cholHostID,cholPartnerID);
                    notAcceptedCholSwaps++;
                }
            }
        }
    }
}





int MCHost::findLipidPairCholNeighbours(int ID1, int ID2)
{   
    std::array<int,4> cholNeigh1=getCholNeighOfLipid(ID1);
    std::array<int,4> cholNeigh2=getCholNeighOfLipid(ID2);
    int numberCholNeigh=0;
    
    for(int i=0;i<4;i++)
    {
        numberCholNeigh+=cholesterinsystem.chols[cholNeigh1[i]].occupied;
        if (!IDinArrayLen4(cholNeigh2[i],cholNeigh1))
        {
            numberCholNeigh+=cholesterinsystem.chols[cholNeigh2[i]].occupied;
        }
    }
//     if (numberCholNeigh==6) std::cout<<"6 lipid pair"<<hostID<<std::endl;
    return numberCholNeigh;
}

int MCHost::findLipidCholPairCholNeighbours(int lipidID, int cholID)
{   
    std::array<int,4> lipidCholNeigh=getCholNeighOfLipid(lipidID);
    std::array<int,4> cholCholNeigh=getCholNeighOfChol(cholID);
    int numberCholNeigh=0;
    
    for(int i=0;i<4;i++)
    {
        numberCholNeigh+=cholesterinsystem.chols[cholCholNeigh[i]].occupied;
        if (!IDinArrayLen4(lipidCholNeigh[i],cholCholNeigh) and lipidCholNeigh[i] != cholID)
        {
            numberCholNeigh+=cholesterinsystem.chols[lipidCholNeigh[i]].occupied;
        }
    }
//     if (numberCholNeigh==5) std::cout<<"5 chol lipid pair"<<hostID <<std::endl;
    return numberCholNeigh;
}


int MCHost::findCholPairCholNeighbours(int ID1, int ID2)
{
    std::array<int,4> cholNeigh1=getCholNeighOfChol(ID1);
    std::array<int,4> cholNeigh2=getCholNeighOfChol(ID2);
    int numberCholNeigh=0;
    
    for(int i=0;i<4;i++)
    {
        if (cholNeigh1[i]!= ID2)
            numberCholNeigh+=cholesterinsystem.chols[cholNeigh1[i]].occupied;
        if (!IDinArrayLen4(cholNeigh2[i],cholNeigh1) and cholNeigh2[i] != ID1)
            numberCholNeigh+=cholesterinsystem.chols[cholNeigh2[i]].occupied;
    }
//     if (numberCholNeigh==6) std::cout<<"6 chol pair"<<hostID<<std::endl;
    return numberCholNeigh;
}



double MCHost::calcSwapEnthalpy()
{
    #ifndef NDEBUG
    std::cout<<"MCHost::calcSwapEnthalpy"<<std::endl;
    #endif

    double H=0;

    std::array<int,4> lipidHostCholNeighbours=getCholNeighOfLipid(hostID);  
    std::array<int,4> lipidPartnerCholNeighbours=getCholNeighOfLipid(partnerID);  

    
    
    for(int i=0;i<4;i++)//H^LC
    {
        if (cholesterinsystem.chols[lipidHostCholNeighbours[i]].occupied)//H^LC host
            H+=lipidproperties->lipidCholEnergieFunction[HOST_LIPID.getType()][findLipidCholPairCholNeighbours(hostID, lipidHostCholNeighbours[i])][HOST_LIPID.getOrder()]*lipidproperties->cholLipidNeigh[getNumberCholNeighOfChol(lipidHostCholNeighbours[i])]/4;
        if (cholesterinsystem.chols[lipidPartnerCholNeighbours[i]].occupied)//H^LC partner
            H+=lipidproperties->lipidCholEnergieFunction[PARTNER_LIPID.getType()][findLipidCholPairCholNeighbours(partnerID, lipidPartnerCholNeighbours[i])][PARTNER_LIPID.getOrder()]*lipidproperties->cholLipidNeigh[getNumberCholNeighOfChol(lipidPartnerCholNeighbours[i])]/4;
    }
    
    
    for(int i=0;i<4;i++)//H^LL
    {
        if (hostNeighbours[i]!=partnerID)
            H+=lipidsystem.calcPairEnthalpy(hostID,hostNeighbours[i],findLipidPairCholNeighbours(hostID,hostNeighbours[i])); //H^LL host     
        if (partnerNeighbours[i]!=hostID)
            H+=lipidsystem.calcPairEnthalpy(partnerID,partnerNeighbours[i],findLipidPairCholNeighbours(partnerID,partnerNeighbours[i])); //H^LL partner
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
    std::cout<<"hostID "<<hostID<<" HOST_posX "<<HOST_posX<<" HOST_posY "<<HOST_posY<<std::endl;
    #endif
    double G=0;
    
    std::array<int,4> lipidHostCholNeighbours=getCholNeighOfLipid(hostID);  

    for(int i=0;i<4;i++)//H^LC
    {
        if (cholesterinsystem.chols[lipidHostCholNeighbours[i]].occupied)//H^LC host
            G+=lipidproperties->lipidCholEnergieFunction[HOST_LIPID.getType()][findLipidCholPairCholNeighbours(hostID, lipidHostCholNeighbours[i])][HOST_LIPID.getOrder()]*lipidproperties->cholLipidNeigh[getNumberCholNeighOfChol(lipidHostCholNeighbours[i])]/4;
    }    
        
    for(int i=0;i<4;i++)
    {
        G+=lipidsystem.calcPairEnthalpy(hostID,hostNeighbours[i],findLipidPairCholNeighbours(hostID,hostNeighbours[i]));
    }
   
    #ifndef NDEBUG
    std::cout<<"H "<<G<<std::endl;
    std::cout<<"kB T S "<<-inputfile->kBT*lipidproperties->entropyFunction[HOST_LIPID.getType()][HOST_LIPID.getOrder()]<<std::endl;
    std::cout<<"self E "<<lipidproperties->selfEnergieFunction[HOST_LIPID.getType()][HOST_LIPID.getOrder()]<<std::endl;
    #endif

    G+=lipidproperties->selfEnergieFunction[HOST_LIPID.getType()][HOST_LIPID.getOrder()]-inputfile->kBT*lipidproperties->entropyFunction[HOST_LIPID.getType()][HOST_LIPID.getOrder()];
    
    
    return G;
}

double MCHost::calcCholSwapEnerg()
{
    #ifndef NDEBUG
    std::cout<<"MCHost::calcCholSwapEnerg"<<std::endl;
    #endif

    double E=0;
    std::array<int,4> cholHostLipidNeighbours=getLipidNeighOfChol(cholHostID);  
    std::array<int,4> cholPartnerLipidNeighbours=getLipidNeighOfChol(cholPartnerID);  
    

    
    for(int i=0;i<4;i++) //change in H^LL next to host
    {
        if (!IDinArrayLen4(cholHostLipidNeighbours[i],cholPartnerLipidNeighbours))
        {
            std::array<int,4> cholHostLipidNeighbourLipidNeighbours=getLipidNeighOfLipid(cholHostLipidNeighbours[i]);
            for(int j=0;j<4;j++)
            {
                if (!IDinArrayLen4(cholHostLipidNeighbourLipidNeighbours[j],cholPartnerLipidNeighbours))
                {
                    if (IDinArrayLen4(cholHostLipidNeighbourLipidNeighbours[j],cholHostLipidNeighbours))
                    {
                        E+=0.5*lipidsystem.calcPairEnthalpy(cholHostLipidNeighbours[i],cholHostLipidNeighbourLipidNeighbours[j],findLipidPairCholNeighbours(cholHostLipidNeighbours[i],cholHostLipidNeighbourLipidNeighbours[j]));
                    }
                    else
                    {
                        E+=lipidsystem.calcPairEnthalpy(cholHostLipidNeighbours[i],cholHostLipidNeighbourLipidNeighbours[j],findLipidPairCholNeighbours(cholHostLipidNeighbours[i],cholHostLipidNeighbourLipidNeighbours[j]));
                    }
                }
            }
        }
    }
    
   
    
    for(int i=0;i<4;i++)//change in H^LL next to partner
    {
        if (IDinArrayLen4(cholPartnerLipidNeighbours[i],cholHostLipidNeighbours))
        {
            std::array<int,4> cholPartnerLipidNeighbourLipidNeighbours=getLipidNeighOfLipid(cholPartnerLipidNeighbours[i]);
            for(int j=0;j<4;j++)
            {
                if (IDinArrayLen4(cholPartnerLipidNeighbourLipidNeighbours[j],cholHostLipidNeighbours))
                {
                    if (IDinArrayLen4(cholPartnerLipidNeighbourLipidNeighbours[j],cholPartnerLipidNeighbours))
                        E+=0.5*lipidsystem.calcPairEnthalpy(cholPartnerLipidNeighbours[i],cholPartnerLipidNeighbourLipidNeighbours[j],findLipidPairCholNeighbours(cholPartnerLipidNeighbours[i],cholPartnerLipidNeighbourLipidNeighbours[j]));
                    else
                        E+=lipidsystem.calcPairEnthalpy(cholPartnerLipidNeighbours[i],cholPartnerLipidNeighbourLipidNeighbours[j],findLipidPairCholNeighbours(cholPartnerLipidNeighbours[i],cholPartnerLipidNeighbourLipidNeighbours[j]));
                }
            }
        }
    }
    
    
    std::array<int,4> cholHostCholNeighbours=getCholNeighOfChol(cholHostID);  
    for(int i=0;i<4;i++)  //H^CC
    {
        if (cholesterinsystem.chols[cholHostCholNeighbours[i]].occupied)
        {
            E+=lipidproperties->cholCholEnergie[findCholPairCholNeighbours(cholHostID,cholHostCholNeighbours[i])]/2;
        }
    }
    
    for(int i=0;i<4;i++) //H^LC
    {
        if (!IDinArrayLen4(cholHostLipidNeighbours[i],cholPartnerLipidNeighbours))
            E+=lipidproperties->lipidCholEnergieFunction[lipidsystem.lipids[cholHostLipidNeighbours[i]].getType()][findLipidCholPairCholNeighbours(cholHostLipidNeighbours[i], cholHostID)][lipidsystem.lipids[cholHostLipidNeighbours[i]].getOrder()]*lipidproperties->cholLipidNeigh[getNumberCholNeighOfChol(cholHostID)]/4;
    }
    

    return E;
}

std::array<int,4> MCHost::getLipidNeighOfChol(int ID)
{
    std::array<int,4> lipidNeighOfChol;
    
    lipidNeighOfChol[0]=lipidsystem.map[cholesterinsystem.chols[ID].posX][cholesterinsystem.chols[ID].posY];
    lipidNeighOfChol[1]=lipidsystem.map[(cholesterinsystem.chols[ID].posX-1+inputfile->width)%inputfile->width][cholesterinsystem.chols[ID].posY];
    lipidNeighOfChol[2]=lipidsystem.map[cholesterinsystem.chols[ID].posX][(cholesterinsystem.chols[ID].posY-1+inputfile->height)%inputfile->height];
    lipidNeighOfChol[3]=lipidsystem.map[(cholesterinsystem.chols[ID].posX-1+inputfile->width)%inputfile->width][(cholesterinsystem.chols[ID].posY-1+inputfile->height)%inputfile->height];
    return lipidNeighOfChol;
}

std::array<int,4> MCHost::getLipidNeighOfLipid(int ID)
{
    std::array<int,4> lipidNeighOfLipid;
    lipidNeighOfLipid[0]=lipidsystem.map[(lipidsystem.lipids[ID].posX+1)%inputfile->width][lipidsystem.lipids[ID].posY];
    lipidNeighOfLipid[1]=lipidsystem.map[(lipidsystem.lipids[ID].posX-1+inputfile->width)%inputfile->width][lipidsystem.lipids[ID].posY];
    lipidNeighOfLipid[2]=lipidsystem.map[lipidsystem.lipids[ID].posX][(lipidsystem.lipids[ID].posY+1)%inputfile->height];
    lipidNeighOfLipid[3]=lipidsystem.map[lipidsystem.lipids[ID].posX][(lipidsystem.lipids[ID].posY-1+inputfile->height)%inputfile->height];
    return lipidNeighOfLipid;
}

std::array<int,4> MCHost::getCholNeighOfLipid(int ID)
{
    std::array<int,4> CholNeighOfLipid;
    CholNeighOfLipid[0]=cholesterinsystem.map[lipidsystem.lipids[ID].posX][lipidsystem.lipids[ID].posY];
    CholNeighOfLipid[1]=cholesterinsystem.map[(lipidsystem.lipids[ID].posX+1)%inputfile->width][lipidsystem.lipids[ID].posY];
    CholNeighOfLipid[2]=cholesterinsystem.map[lipidsystem.lipids[ID].posX][(lipidsystem.lipids[ID].posY+1)%inputfile->height];
    CholNeighOfLipid[3]=cholesterinsystem.map[(lipidsystem.lipids[ID].posX+1)%inputfile->width][(lipidsystem.lipids[ID].posY+1)%inputfile->height];
    return CholNeighOfLipid;
}

std::array<int,4> MCHost::getCholNeighOfChol(int ID)
{
    std::array<int,4> CholNeighOfChol;
    CholNeighOfChol[0]=cholesterinsystem.map[(cholesterinsystem.chols[ID].posX+1)%inputfile->width][cholesterinsystem.chols[ID].posY];
    CholNeighOfChol[1]=cholesterinsystem.map[(cholesterinsystem.chols[ID].posX-1+inputfile->width)%inputfile->width][cholesterinsystem.chols[ID].posY];
    CholNeighOfChol[2]=cholesterinsystem.map[cholesterinsystem.chols[ID].posX][(cholesterinsystem.chols[ID].posY+1)%inputfile->height];
    CholNeighOfChol[3]=cholesterinsystem.map[cholesterinsystem.chols[ID].posX][(cholesterinsystem.chols[ID].posY-1+inputfile->height)%inputfile->height];
    return CholNeighOfChol;
}

int MCHost::getNumberCholNeighOfLipid(int ID)
{
    int NumberCholNeighOfLipid=0;   
    NumberCholNeighOfLipid+=cholesterinsystem.chols[cholesterinsystem.map[lipidsystem.lipids[ID].posX][lipidsystem.lipids[ID].posY]].occupied;
    NumberCholNeighOfLipid+=cholesterinsystem.chols[cholesterinsystem.map[(lipidsystem.lipids[ID].posX+1)%inputfile->width][lipidsystem.lipids[ID].posY]].occupied;
    NumberCholNeighOfLipid+=cholesterinsystem.chols[cholesterinsystem.map[lipidsystem.lipids[ID].posX][(lipidsystem.lipids[ID].posY+1)%inputfile->height]].occupied;
    NumberCholNeighOfLipid+=cholesterinsystem.chols[cholesterinsystem.map[(lipidsystem.lipids[ID].posX+1)%inputfile->width][(lipidsystem.lipids[ID].posY+1)%inputfile->height]].occupied;
    return NumberCholNeighOfLipid;
}

int MCHost::getNumberCholNeighOfChol(int ID)
{
    int NumberCholNeighOfChol=0;
    NumberCholNeighOfChol+=cholesterinsystem.chols[cholesterinsystem.map[(cholesterinsystem.chols[ID].posX+1)%inputfile->width][cholesterinsystem.chols[ID].posY]].occupied;
    NumberCholNeighOfChol+=cholesterinsystem.chols[cholesterinsystem.map[(cholesterinsystem.chols[ID].posX-1+inputfile->width)%inputfile->width][cholesterinsystem.chols[ID].posY]].occupied;
    NumberCholNeighOfChol+=cholesterinsystem.chols[cholesterinsystem.map[cholesterinsystem.chols[ID].posX][(cholesterinsystem.chols[ID].posY+1)%inputfile->height]].occupied;
    NumberCholNeighOfChol+=cholesterinsystem.chols[cholesterinsystem.map[cholesterinsystem.chols[ID].posX][(cholesterinsystem.chols[ID].posY-1+inputfile->height)%inputfile->height]].occupied;
    return NumberCholNeighOfChol;
}


bool MCHost::IDinArrayLen4(int& ID, std::array<int,4>& Array )
{
    for(int i=0;i<4;i++) 
        if (ID==Array[i])
            return true;
    return false;
}
void MCHost::setHost(int x, int y)
{
    hostID=lipidsystem.map[x][y];
    hostNeighbours=getLipidNeighOfLipid(hostID);
}

void MCHost::setHost(int ID)
{
    hostID=ID;
    hostNeighbours=getLipidNeighOfLipid(hostID);
}

void MCHost::setRNDHost()
{
    hostID=enhance::random_int(0,inputfile->height*inputfile->width-1);
    hostNeighbours=getLipidNeighOfLipid(hostID);
}

void MCHost::setPartner()
{
    partnerID=hostNeighbours[enhance::random_int(0,3)];
    partnerNeighbours=getLipidNeighOfLipid(partnerID);
}


bool MCHost::setRNDCholHost()
{
    cholHostID=enhance::random_int(0,inputfile->height*inputfile->width-1);
    return cholesterinsystem.chols[cholHostID].occupied;
}

bool MCHost::setCholHost(int ID)
{
    cholHostID=ID;
    return cholesterinsystem.chols[cholHostID].occupied;
}

bool MCHost::setCholPartner()
{
    std::array<int,4> cholHostCholNeighbours=getCholNeighOfChol(cholHostID);
    cholPartnerID=cholHostCholNeighbours[enhance::random_int(0,3)];
    return cholesterinsystem.chols[cholPartnerID].occupied;
}


bool MCHost::acceptance(const double FreeEnergie1, const double FreeEnergie2)
{
    #ifndef NDEBUG
    std::cout<<"MCHost::acceptance:   DeltaG: "<<(FreeEnergie2-FreeEnergie1)<<" DeltaG/kbT: " <<(FreeEnergie2-FreeEnergie1)/inputfile->kBT<<std::endl;
    #endif

    return enhance::random_double(0.0, 1.0) < enhance::fastExp((FreeEnergie1-FreeEnergie2)/inputfile->kBT);
    
}
