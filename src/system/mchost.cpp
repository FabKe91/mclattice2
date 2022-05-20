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
    std::cerr<<"MCHost::setup"<<std::endl;
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
    std::cerr<<"MCHost::setupForRestart"<<std::endl;
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
    std::cerr<<"MCHost::run"<<std::endl;
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
    std::cerr<<"MCHost::doSystemloop"<<std::endl;
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




std::array<int, 3> MCHost::findLipidPairCholNeighbours(int ID1, int ID2)
{   
    std::array<int,4> cholNeigh1=getCholNeighOfLipid(ID1);
    std::array<int,4> cholNeigh2=getCholNeighOfLipid(ID2);
    //int numberCholNeigh=0;

    int mutualNc = 0;
    int hostNc = 0;
    int neibNc = 0;
    
    for(int i=0;i<4;i++)
    {
        hostNc+=cholesterinsystem.chols[cholNeigh1[i]].occupied;
        neibNc+=cholesterinsystem.chols[cholNeigh2[i]].occupied;

        if (IDinArrayLen4(cholNeigh2[i],cholNeigh1))
        {
            mutualNc+=cholesterinsystem.chols[cholNeigh2[i]].occupied;
        }
    }

    //     if (numberCholNeigh==6) std::cout<<"6 lipid pair"<<hostID<<std::endl;

    LipidPairCholNeighbors[0] = hostNc;
    LipidPairCholNeighbors[1] = neibNc;
    LipidPairCholNeighbors[2] = mutualNc;

    #ifndef NDEBUG
    std::cerr<<"MCHost::findLipidPairCholNeighbours"<<std::endl;
    std::cerr<<"hostNc: "<<hostNc<<" neibNc: "<<neibNc<<" mutualNc: "<<mutualNc<<std::endl;
    #endif 

    return LipidPairCholNeighbors;
}

std::array<int, 3> MCHost::findLipidCholPairCholNeighbours(int lipidID, int cholID)
{   
    std::array<int,4> lipidCholNeigh=getCholNeighOfLipid(lipidID);
    std::array<int,4> cholCholNeigh=getCholNeighOfChol(cholID);
    //int numberCholNeigh=0;

    int mutualNc = 0;
    int hostNc = 0;
    int neibNc = 0;
    
    for(int i=0;i<4;i++)
    {
        neibNc+=cholesterinsystem.chols[cholCholNeigh[i]].occupied;

        if ( lipidCholNeigh[i] != cholID )
        {
            hostNc+=cholesterinsystem.chols[lipidCholNeigh[i]].occupied;
        }

        if (IDinArrayLen4(lipidCholNeigh[i],cholCholNeigh) and lipidCholNeigh[i] != cholID)
        {
            //numberCholNeigh+=cholesterinsystem.chols[lipidCholNeigh[i]].occupied;
            mutualNc+=cholesterinsystem.chols[lipidCholNeigh[i]].occupied;
        }

    }

    LipidCholPairCholNeighbors[0] = hostNc;
    LipidCholPairCholNeighbors[1] = neibNc;
    LipidCholPairCholNeighbors[2] = mutualNc;

    #ifndef NDEBUG
    std::cerr<<"MCHost::findLipidCholPairCholNeighbours"<<std::endl;
    std::cerr<<" cholID "<<cholID<<" hostNc: "<<hostNc<<" neibNc: "<<neibNc<<" mutualNc "<<mutualNc<<std::endl;
    for(int i=0;i<4;i++) {std::cerr<<"cc-"<<i<<":"<<cholCholNeigh[i]<<"="<<cholesterinsystem.chols[cholCholNeigh[i]].occupied<<" ";}
    std::cerr<<std::endl;
    for(int i=0;i<4;i++) {std::cerr<<"lc-"<<i<<":"<<lipidCholNeigh[i]<<"="<<cholesterinsystem.chols[lipidCholNeigh[i]].occupied<<" ";}
    std::cerr<<std::endl;
    #endif 

    //if (numberCholNeigh==5) std::cout<<"5 chol lipid pair"<<hostID <<std::endl;
    return LipidCholPairCholNeighbors;
}


std::array<int, 3> MCHost::findCholPairCholNeighbours(int ID1, int ID2)
{
    std::array<int,4> cholNeigh1=getCholNeighOfChol(ID1);
    std::array<int,4> cholNeigh2=getCholNeighOfChol(ID2);
    //int numberCholNeigh=0;

    int mutualNc = 0; // CholPairs cannot have mutual Nc
    int hostNc = 0;
    int neibNc = 0;
    
    for(int i=0;i<4;i++)
    {
        if (cholNeigh1[i] != ID2)
            hostNc+=cholesterinsystem.chols[cholNeigh1[i]].occupied;
        if (cholNeigh2[i] != ID1)
            neibNc+=cholesterinsystem.chols[cholNeigh2[i]].occupied;
    }

    CholPairCholNeighbors[0] = hostNc;
    CholPairCholNeighbors[1] = neibNc;
    CholPairCholNeighbors[2] = mutualNc;

    #ifndef NDEBUG
    std::cerr<<"MCHost::findCholPairCholNeighbours"<<std::endl;
    std::cerr<<"hostNc:"<<hostNc<<"neibNc:"<<neibNc<<"mutualNc"<<mutualNc<<std::endl;
    #endif 
    //if (numberCholNeigh==6) std::cout<<"6 chol pair"<<hostID<<std::endl;
    return CholPairCholNeighbors;
}

double MCHost::calcSwapEnthalpy()
{
    #ifndef NDEBUG
    std::cerr<<"MCHost::calcSwapEnthalpy"<<std::endl;
    #endif

    double H=0;

    std::array<int,4> lipidHostCholNeighbours=getCholNeighOfLipid(hostID);  
    std::array<int,4> lipidPartnerCholNeighbours=getCholNeighOfLipid(partnerID);  

    
    for(int i=0;i<4;i++)//H^LC
    {
        std::array<int, 3> hostCholPairCholNeighbors=findLipidCholPairCholNeighbours(hostID, lipidHostCholNeighbours[i]);
        std::array<int, 3> partnerCholPairCholNeighbors=findLipidCholPairCholNeighbours(partnerID, lipidPartnerCholNeighbours[i]);
        int hostCholPairCholNeighbors_tot = hostCholPairCholNeighbors[0] + hostCholPairCholNeighbors[1] - hostCholPairCholNeighbors[2];
        int partnerCholPairCholNeighbors_tot = partnerCholPairCholNeighbors[0] + partnerCholPairCholNeighbors[1] - partnerCholPairCholNeighbors[2];

        if (cholesterinsystem.chols[lipidHostCholNeighbours[i]].occupied)//H^LC host
            // H += LC[h_i, c_i](Nc, S) * N[h_i, c_i](Nc, T) / 4
            H+=lipidproperties->lipidCholEnergieFunction[HOST_LIPID.getType()][hostCholPairCholNeighbors_tot][HOST_LIPID.getOrder()] * lipidproperties->cholLipidNeigh[HOST_LIPID.getType()][getNumberCholNeighOfChol(lipidHostCholNeighbours[i])] / 4;

        if (cholesterinsystem.chols[lipidPartnerCholNeighbours[i]].occupied)//H^LC partner
            H+=lipidproperties->lipidCholEnergieFunction[PARTNER_LIPID.getType()][partnerCholPairCholNeighbors_tot][PARTNER_LIPID.getOrder()] * lipidproperties->cholLipidNeigh[PARTNER_LIPID.getType()][getNumberCholNeighOfChol(lipidPartnerCholNeighbours[i])] / 4;
    }
    
    
    for(int i=0;i<4;i++)//H^LL
    {
        if (hostNeighbours[i]!=partnerID)
            H+=lipidsystem.calcPairEnthalpy(hostID,hostNeighbours[i],findLipidPairCholNeighbours(hostID,hostNeighbours[i])); //H^LL host     
        if (partnerNeighbours[i]!=hostID)
            H+=lipidsystem.calcPairEnthalpy(partnerID,partnerNeighbours[i],findLipidPairCholNeighbours(partnerID,partnerNeighbours[i])); //H^LL partner
    }   

    #ifndef NDEBUG
    std::cerr<<"H "<<H<<std::endl;
    #endif
    
    return H;
}






double MCHost::calcHostFreeEnerg()
{
    #ifndef NDEBUG
    std::cerr<<"MCHost::calcHostFreeEnerg"<<std::endl;
    std::cerr<<"hostID "<<hostID<<" HOST_posX "<<HOST_posX<<" HOST_posY "<<HOST_posY<<std::endl;
    #endif
    double G=0;
    
    std::array<int,4> lipidHostCholNeighbours=getCholNeighOfLipid(hostID);  

    for(int i=0;i<4;i++)//H^LC
    {

        std::array<int, 3> LipidCholPairCholNeighbors=findLipidCholPairCholNeighbours(hostID, lipidHostCholNeighbours[i]);
        int LipidCholPairCholNeighbors_tot = LipidCholPairCholNeighbors[0] + LipidCholPairCholNeighbors[1] - LipidCholPairCholNeighbors[2];

        #ifndef NDEBUG
        std::cerr<<"LipidCholPairCholNeighbors_tot"<<LipidCholPairCholNeighbors_tot<<std::endl;
        #endif

        if (cholesterinsystem.chols[lipidHostCholNeighbours[i]].occupied)//H^LC host
            G+=lipidproperties->lipidCholEnergieFunction[HOST_LIPID.getType()][LipidCholPairCholNeighbors_tot][HOST_LIPID.getOrder()]*lipidproperties->cholLipidNeigh[HOST_LIPID.getType()][getNumberCholNeighOfChol(lipidHostCholNeighbours[i])]/4;
    }    
        
    for(int i=0;i<4;i++)//H^LL
    {
        G+=lipidsystem.calcPairEnthalpy(hostID,hostNeighbours[i],findLipidPairCholNeighbours(hostID,hostNeighbours[i]));
    }
   
    #ifndef NDEBUG
    std::cerr<<"H "<<G<<std::endl;
    std::cerr<<"kB T S "<<-inputfile->kBT*lipidproperties->entropyFunction[HOST_LIPID.getType()][HOST_LIPID.getOrder()]<<std::endl;
    std::cerr<<"self E "<<lipidproperties->selfEnergieFunction[HOST_LIPID.getType()][HOST_LIPID.getOrder()]<<std::endl;
    #endif

    G+=lipidproperties->selfEnergieFunction[HOST_LIPID.getType()][HOST_LIPID.getOrder()]-inputfile->kBT*lipidproperties->entropyFunction[HOST_LIPID.getType()][HOST_LIPID.getOrder()];
    
    
    return G;
}

double MCHost::calcCholSwapEnerg()
{
    #ifndef NDEBUG
    std::cerr<<"MCHost::calcCholSwapEnerg"<<std::endl;
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
        std::array<int, 3> cholpairCholNeighbors = findCholPairCholNeighbours(cholHostID,cholHostCholNeighbours[i]);
        int cholpairCholNeighbors_tot = cholpairCholNeighbors[0] + cholpairCholNeighbors[1] - cholpairCholNeighbors[2];
        if (cholesterinsystem.chols[cholHostCholNeighbours[i]].occupied)
        {
            E+=lipidproperties->cholCholEnergie[cholpairCholNeighbors_tot];
        }
    }
    
    for(int i=0;i<4;i++) //H^LC
    {
        std::array<int, 3> lipidCholPairCholNeighbours = findLipidCholPairCholNeighbours(cholHostLipidNeighbours[i], cholHostID);
        int lipidCholPairCholNeighbours_tot = lipidCholPairCholNeighbours[0] + lipidCholPairCholNeighbours[1] - lipidCholPairCholNeighbours[2];
        if (!IDinArrayLen4(cholHostLipidNeighbours[i],cholPartnerLipidNeighbours))
            E+=lipidproperties->lipidCholEnergieFunction[lipidsystem.lipids[cholHostLipidNeighbours[i]].getType()][lipidCholPairCholNeighbours_tot][lipidsystem.lipids[cholHostLipidNeighbours[i]].getOrder()] * lipidproperties->cholLipidNeigh[lipidsystem.lipids[cholHostLipidNeighbours[i]].getType()][getNumberCholNeighOfChol(cholHostID)] / 4;
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

//int MCHost::getNumberCholNeighOfLipid(int ID)
//{
//    int NumberCholNeighOfLipid=0;   
//    NumberCholNeighOfLipid+=cholesterinsystem.chols[cholesterinsystem.map[lipidsystem.lipids[ID].posX][lipidsystem.lipids[ID].posY]].occupied;
//    NumberCholNeighOfLipid+=cholesterinsystem.chols[cholesterinsystem.map[(lipidsystem.lipids[ID].posX+1)%inputfile->width][lipidsystem.lipids[ID].posY]].occupied;
//    NumberCholNeighOfLipid+=cholesterinsystem.chols[cholesterinsystem.map[lipidsystem.lipids[ID].posX][(lipidsystem.lipids[ID].posY+1)%inputfile->height]].occupied;
//    NumberCholNeighOfLipid+=cholesterinsystem.chols[cholesterinsystem.map[(lipidsystem.lipids[ID].posX+1)%inputfile->width][(lipidsystem.lipids[ID].posY+1)%inputfile->height]].occupied;
//    return NumberCholNeighOfLipid;
//}

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
    std::cerr<<"MCHost::acceptance:   DeltaG: "<<(FreeEnergie2-FreeEnergie1)<<" DeltaG/kbT: " <<(FreeEnergie2-FreeEnergie1)/inputfile->kBT<<std::endl;
    #endif

    return enhance::random_double(0.0, 1.0) < enhance::fastExp((FreeEnergie1-FreeEnergie2)/inputfile->kBT);
    
}
