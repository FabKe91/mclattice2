#include "lipidsystem.h"

Lipidsystem::Lipidsystem()
{

}

void Lipidsystem::readParas(std::shared_ptr<InputFile> _inputfile,    std::shared_ptr<LipidProperties> _lipidproperties)
{
    lipidproperties=_lipidproperties;
    inputfile=_inputfile;
    height=inputfile->paras["height"];
    width=inputfile->paras["width"];
    kBT=inputfile->paras["kB"]*inputfile->paras["T"];

    
}


void Lipidsystem::setup()
{   
    #ifndef NDEBUG
    std::cout<<"Lipidsystem::setup"<<std::endl;
    #endif
    
    
    lastSwappedPos = new unsigned int*[2];
    lastSwappedPos[0]= new unsigned int[2];
    lastSwappedPos[1]= new unsigned int[2];
    
    map =new unsigned int*[width];
    for(int i=0;i<width;i++)
    {
        map[i] =new unsigned int[height];
        for(int j=0;j<height;j++)
        {   
            map[i][j]=i*height+j;
            lipids.push_back(Lipid(i*height+j,0,0));
//             lipids.push_back(Lipid(i*10+j,0,enhance::random_int(0,inputfile->types[0].maxOrder)));
        }
    }
}


int Lipidsystem::getMeanOrder()
{
    #ifndef NDEBUG
    std::cout<<"Lipidsystem::getMeanOrder"<<std::endl;
    #endif
    
    int mean=0;
    int count=0;
    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
        {   
            count++;
            mean += lipids[map[i][j]].getOrderPara();
        }
    }
    return mean/count;
}

std::vector<int> Lipidsystem::getOrderDestr()
{
    std::vector<int> destr((int)inputfile->paras["maxOrderIndex"]+1,0);
    int count=0;
    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
        {   
            count++;
            destr[lipids[map[i][j]].getOrderPara()]++;
        }
    }
    
}


const boost::multi_array<int,2> Lipidsystem::getOrderParas()
{
//     std::cout<<"getOrderParas"<<width<<height<<std::endl;

    boost::multi_array<int,2> data(boost::extents[width][height]);
    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
        {   
                data[i][j] = lipids[map[i][j]].getOrderPara();
        }
    }
    return data;
}

const boost::multi_array<int,2> Lipidsystem::getTypes()
{
//     std::cout<<"getTypes"<<width<<height<<std::endl;

    boost::multi_array<int,2> data(boost::extents[width][height]);
    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
        {   
                data[i][j] = lipids[map[i][j]].getType();
        }
    }
    return data;
}



void Lipidsystem::setHost(int i, int j)
{
    lastSwappedPos[0][0]=i;
    lastSwappedPos[0][1]=j;
    lastSwappedIDs[0]=map[i][j];
}

void Lipidsystem::setRNDHost()
{
    lastSwappedPos[0][0]=enhance::random_uns_int(0,width-1);
    lastSwappedPos[0][1]=enhance::random_uns_int(0,height-1);
    lastSwappedIDs[0]=map[lastSwappedPos[0][0]][lastSwappedPos[0][1]];
}

void Lipidsystem::setPartner()
{
    rdnPartnerNumber=enhance::random_int(0,3);
    
    if (rdnPartnerNumber==0)
    {
        lastSwappedPos[1][0]=(lastSwappedPos[0][0]+1)%width;
        lastSwappedPos[1][1]=lastSwappedPos[0][0];
    }    
    else if (rdnPartnerNumber==1)
    {
        lastSwappedPos[1][0]=(lastSwappedPos[0][0]-1)%width;
        lastSwappedPos[1][1]=lastSwappedPos[0][0];
    }    
    else if (rdnPartnerNumber==2)
    {
        lastSwappedPos[1][0]=lastSwappedPos[0][0];
        lastSwappedPos[1][1]=(lastSwappedPos[0][0]+1)%width;
    }
    else
    {
        lastSwappedPos[1][0]=lastSwappedPos[0][0];
        lastSwappedPos[1][1]=(lastSwappedPos[0][0]-1)%width;
    }
    
    lastSwappedIDs[1]=map[lastSwappedPos[1][0]][lastSwappedPos[1][1]];
}

void Lipidsystem::swap()
{
    #ifndef NDEBUG
    std::cout<<"Lipidsystem::swap"<<std::endl;
    #endif

    map[lastSwappedPos[0][0]][lastSwappedPos[0][1]]=lastSwappedIDs[1];
    map[lastSwappedPos[1][0]][lastSwappedPos[1][1]]=lastSwappedIDs[0];
    int IDbuffer=lastSwappedIDs[1];
    lastSwappedIDs[1]=lastSwappedIDs[0];
    lastSwappedIDs[0]=IDbuffer;
}



double Lipidsystem::calcSwapEnthalpy()
{
    #ifndef NDEBUG
    std::cout<<"Lipidsystem::calcSwapEnthalpy"<<std::endl;
    #endif

    double H=0;
    
    #define ID0 lastSwappedIDs[0]
    #define posX0 lastSwappedPos[0][0]
    #define posY0 lastSwappedPos[0][1]
    #define ID1 lastSwappedIDs[1]
    #define posX1 lastSwappedPos[1][0]
    #define posY1 lastSwappedPos[1][1]


    if (rdnPartnerNumber==0)
    {
        H+=calcPairEnthalpy(ID0,map[(posX0-1)%width][posY0]);
        H+=calcPairEnthalpy(ID0,map[posX0][(posY0+1)%height]);
        H+=calcPairEnthalpy(ID0,map[posX0][(posY0-1)%height]);
        
        H+=calcPairEnthalpy(ID1,map[(posX1+1)%width][posY1]);
        H+=calcPairEnthalpy(ID1,map[posX1][(posY1+1)%height]);
        H+=calcPairEnthalpy(ID1,map[posX1][(posY1-1)%height]);
    }    
    else if (rdnPartnerNumber==1)
    {
        H+=calcPairEnthalpy(ID0,map[(posX0+1)%width][posY0]);
        H+=calcPairEnthalpy(ID0,map[posX0][(posY0+1)%height]);
        H+=calcPairEnthalpy(ID0,map[posX0][(posY0-1)%height]);
        
        H+=calcPairEnthalpy(ID1,map[(posX1-1)%width][posY1]);
        H+=calcPairEnthalpy(ID1,map[posX1][(posY1+1)%height]);
        H+=calcPairEnthalpy(ID1,map[posX1][(posY1-1)%height]);
    }    
    else if (rdnPartnerNumber==2)
    {
        H+=calcPairEnthalpy(ID0,map[(posX0+1)%width][posY0]);
        H+=calcPairEnthalpy(ID0,map[(posX0-1)%width][posY0]);
        H+=calcPairEnthalpy(ID0,map[posX0][(posY0-1)%height]);
        
        H+=calcPairEnthalpy(ID1,map[(posX1+1)%width][posY1]);
        H+=calcPairEnthalpy(ID1,map[(posX1-1)%width][posY1]);
        H+=calcPairEnthalpy(ID1,map[posX1][(posY1+1)%height]);
    }
    else
    {
        H+=calcPairEnthalpy(ID0,map[(posX0+1)%width][posY0]);
        H+=calcPairEnthalpy(ID0,map[(posX0-1)%width][posY0]);
        H+=calcPairEnthalpy(ID0,map[posX0][(posY0+1)%height]);
        
        H+=calcPairEnthalpy(ID1,map[(posX1+1)%width][posY1]);
        H+=calcPairEnthalpy(ID1,map[(posX1-1)%width][posY1]);
        H+=calcPairEnthalpy(ID1,map[posX1][(posY1-1)%height]);
    }
    #ifndef NDEBUG
    std::cout<<"H "<<H<<std::endl;
    #endif
    
    return H;
}

double Lipidsystem::calcHostFreeEnerg()
{
    #ifndef NDEBUG
    std::cout<<"Lipidsystem::calcHostFreeEnerg"<<std::endl;
    #endif

    
    #define ID0 lastSwappedIDs[0]
    #define posX0 lastSwappedPos[0][0]
    #define posY0 lastSwappedPos[0][1]

    
    double G=0;

    G+=calcPairEnthalpy(ID0,map[(posX0-1)%width][posY0]);
    G+=calcPairEnthalpy(ID0,map[(posX0+1)%width][posY0]);
    G+=calcPairEnthalpy(ID0,map[posX0][(posY0-1)%height]);
    G+=calcPairEnthalpy(ID0,map[posX0][(posY0-1)%height]);

   
    #ifndef NDEBUG
    std::cout<<"H "<<G<<std::endl;
    std::cout<<"kB T S "<<-kBT*lipidproperties->entropyFunction[lipids[ID0].getType()][lipids[ID0].getOrderPara()]<<std::endl;
    
    #endif

    G-=kBT*lipidproperties->entropyFunction[lipids[ID0].getType()][lipids[ID0].getOrderPara()];
    
    
    return G;
}




double Lipidsystem::calcPairEnthalpy(unsigned int ID_1,unsigned int ID_2)
{ 
    int type1=lipids[ID_1].getType();
    int type2=lipids[ID_2].getType();
    int order1=lipids[ID_1].getOrderPara();
    int order2=lipids[ID_2].getOrderPara();
    
//     #define type1 lipids[ID1].getType()
//     #define type2 lipids[ID2].getType()
//     #define order1 lipids[ID1].getOrderPara()
//     #define order2 lipids[ID2].getOrderPara()
//     
    
    if (type1>=type2)
    {
        return lipidproperties->enthalpyFunction[type1][type2][(int)((order1+order2)/2)]*(lipidproperties->neighbourFunction[type1][order1]+lipidproperties->neighbourFunction[type2][order2])/16;
    }
    else
    {
        return lipidproperties->enthalpyFunction[type2][type1][(int)((order1+order2)/2)]*(lipidproperties->neighbourFunction[type1][order1]+lipidproperties->neighbourFunction[type2][order2])/16;;
    }

}


void Lipidsystem::fluctuate()
{
    #ifndef NDEBUG
    std::cout<<"Lipidsystem::fluctuate"<<std::endl;
    #endif

    oldOrder=lipids[lastSwappedIDs[0]].getOrderPara();
    int maxOrder=inputfile->types[lipids[lastSwappedIDs[0]].getType()].maxOrder;
    int minOrder=inputfile->types[lipids[lastSwappedIDs[0]].getType()].minOrder;
    int maxFluc=inputfile->types[lipids[lastSwappedIDs[0]].getType()].maxFluc;
    
    lipids[lastSwappedIDs[0]].setOrderPara(enhance::random_int((oldOrder-maxFluc+minOrder+std::abs(oldOrder-maxFluc-minOrder))/2,(oldOrder+maxFluc+maxOrder-std::abs(oldOrder+maxFluc-maxOrder))/2));
    
    #ifndef NDEBUG
    std::cout<<"oldOrder: "<<oldOrder<<" maxOrder: "<<maxOrder<<" minOrder: "<<minOrder<<" maxFluc: "<<maxFluc<<std::endl;
    std::cout<<"searching order between "<<(oldOrder-maxFluc+minOrder+std::abs(oldOrder-maxFluc-minOrder))/2<<" "<<(oldOrder+maxFluc+maxOrder-std::abs(oldOrder+maxFluc-maxOrder))/2<<" new order: "<<lipids[lastSwappedIDs[0]].getOrderPara()<<std::endl;
    #endif
    
    
    
}

void Lipidsystem::fluctuateBack()
{
    #ifndef NDEBUG
    std::cout<<"Lipidsystem::fluctuateBack"<<std::endl;
    #endif
    lipids[lastSwappedIDs[0]].setOrderPara(oldOrder);
    
}







