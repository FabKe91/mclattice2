#include "lipidsystem.h"


#define ID0 lastSwappedIDs[0]
#define posX0 lipids[ID0].posX
#define posY0 lipids[ID0].posY
#define ID1 lastSwappedIDs[1]
#define posX1 lipids[ID1].posX
#define posY1 lipids[ID1].posY
   

Lipidsystem::Lipidsystem()
{

}

void Lipidsystem::readParas( std::shared_ptr<LipidProperties> _lipidproperties,std::shared_ptr<InputFile> _inputfile)
{
    lipidproperties=_lipidproperties;
    inputfile=_inputfile;
    height=inputfile->height;
    width=inputfile->width;
}


void Lipidsystem::setup()
{   
    #ifndef NDEBUG
    std::cerr<<"Lipidsystem::setup"<<std::endl;
    #endif
    
    
    int defaultOrderIndex=(inputfile->paras.at("defaultOrderPara")-inputfile->paras.at("minOrder"))/inputfile->paras.at("DeltaOrder");
    
    //setting up type distr
    int N=width*height;
    std::vector<int> types{};
    for(int type=0;type<inputfile->nType-1;type++)
        for(int j=0;j<inputfile->concentrations[type]*N and types.size()<=N ;j++)
            types.push_back(type);

    while (types.size()<N) types.push_back(inputfile->nType-1); 
    std::shuffle(types.begin(), types.end(), enhance::rand_engine); //randomize order of types


    map =new  int*[width];
    for(int i=0;i<width;i++)
    {
        map[i] =new  int[height];
        for(int j=0;j<height;j++)
        {   
            map[i][j]=i*height+j;
            lipids.push_back(Lipid(i*height+j,types[i*height+j],defaultOrderIndex,i,j));
        }
    }
}


int Lipidsystem::getMeanOrder()
{
    #ifndef NDEBUG
    std::cerr<<"Lipidsystem::getMeanOrder"<<std::endl;
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

std::vector<int> Lipidsystem::getOrderDistr()
{
    std::vector<int> distr((int)inputfile->paras.at("maxOrderIndex")+1,0);
    for(int i=0;i<width;i++)
    for(int j=0;j<height;j++)
        distr[lipids[map[i][j]].getOrderPara()]++;
    return distr;
    
}


const boost::multi_array<int,2> Lipidsystem::getOrderParas()
{
//     std::cout<<"getOrderParas"<<width<<height<<std::endl;

    boost::multi_array<int,2> data(boost::extents[width][height]);
    for(int i=0;i<width;i++)    
    for(int j=0;j<height;j++)
        data[i][j] = lipids[map[i][j]].getOrderPara();
    return data;
}

const boost::multi_array<int,2> Lipidsystem::getTypes()
{
//     std::cout<<"getTypes"<<width<<height<<std::endl;

    boost::multi_array<int,2> data(boost::extents[width][height]);
    for(int i=0;i<width;i++)
    for(int j=0;j<height;j++)  
        data[i][j] = lipids[map[i][j]].getType();
    return data;
}

const boost::multi_array<int,2> Lipidsystem::getIDs()
{
//     std::cout<<"getTypes"<<width<<height<<std::endl;

    boost::multi_array<int,2> data(boost::extents[width][height]);
    for(int i=0;i<width;i++)
    for(int j=0;j<height;j++)  
        data[i][j] = lipids[map[i][j]].getID();
    return data;
}



void Lipidsystem::setHost(int x, int y)
{
    ID0=map[x][y];
}

void Lipidsystem::setHost(int ID)
{
    ID0=ID;
}

void Lipidsystem::setRNDHost()
{
    ID0=enhance::random_int(0,height*width);
}


void Lipidsystem::printMap()
{
    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
            std::cout<<map[i][j]<<" ";
        std::cout<<std::endl;
    }
    
}

void Lipidsystem::setPartner()
{
    rdnPartnerNumber=enhance::random_int(0,3);
    
    switch (rdnPartnerNumber)
    {
        case 0 :    ID1=map[(posX0+1)%width][posY0]; //up
                    break;
        case 1 :    ID1=map[(posX0-1+width)%width][posY0]; //down
                    break;
        case 2 :    ID1=map[posX0][(posY0+1)%height]; //right
                    break;
        case 3 :    ID1=map[posX0][(posY0-1+height)%height]; //left
    }
}

void Lipidsystem::swap()
{
    #ifndef NDEBUG
    std::cerr<<"Lipidsystem::swap"<<std::endl;
    #endif
//     std::cout<<posX0<<" "<<posY0<<std::endl;
//     std::cout<<posX1<<" "<<posY1<<std::endl;
    
//         std::cout<<"Lipidsystem::swap"<<std::endl;

    std::swap(posX0,posX1);
    std::swap(posY0,posY1);
    std::swap(map[posX0][posY0],map[posX1][posY1]);
    
//     std::cout<<posX0<<" "<<posY0<<std::endl;
//     std::cout<<posX1<<" "<<posY1<<std::endl;
}



double Lipidsystem::calcSwapEnthalpy()
{
    #ifndef NDEBUG
    std::cerr<<"Lipidsystem::calcSwapEnthalpy"<<std::endl;
    #endif

    double H=0;
    


    switch(rdnPartnerNumber)
    {
        case 0: H+=calcPairEnthalpy(ID0,map[(posX0-1+width)%width][posY0]);
                H+=calcPairEnthalpy(ID0,map[posX0][(posY0+1)%height]);
                H+=calcPairEnthalpy(ID0,map[posX0][(posY0-1+height)%height]);
                
                H+=calcPairEnthalpy(ID1,map[(posX1+1)%width][posY1]);
                H+=calcPairEnthalpy(ID1,map[posX1][(posY1+1)%height]);
                H+=calcPairEnthalpy(ID1,map[posX1][(posY1-1+height)%height]);
                break;
                
        case 1: H+=calcPairEnthalpy(ID0,map[(posX0+1)%width][posY0]);
                H+=calcPairEnthalpy(ID0,map[posX0][(posY0+1)%height]);
                H+=calcPairEnthalpy(ID0,map[posX0][(posY0-1+height)%height]);
                
                H+=calcPairEnthalpy(ID1,map[(posX1-1+width)%width][posY1]);
                H+=calcPairEnthalpy(ID1,map[posX1][(posY1+1)%height]);
                H+=calcPairEnthalpy(ID1,map[posX1][(posY1-1+height)%height]);
                break;
       
        case 2: H+=calcPairEnthalpy(ID0,map[(posX0+1)%width][posY0]);
                H+=calcPairEnthalpy(ID0,map[(posX0-1+width)%width][posY0]);
                H+=calcPairEnthalpy(ID0,map[posX0][(posY0-1+height)%height]);
                
                H+=calcPairEnthalpy(ID1,map[(posX1+1)%width][posY1]);
                H+=calcPairEnthalpy(ID1,map[(posX1-1)%width][posY1]);
                H+=calcPairEnthalpy(ID1,map[posX1][(posY1+1)%height]);
                break;
                
        case 3: H+=calcPairEnthalpy(ID0,map[(posX0+1)%width][posY0]);
                H+=calcPairEnthalpy(ID0,map[(posX0-1+width)%width][posY0]);
                H+=calcPairEnthalpy(ID0,map[posX0][(posY0+1)%height]);
                
                H+=calcPairEnthalpy(ID1,map[(posX1+1)%width][posY1]);
                H+=calcPairEnthalpy(ID1,map[(posX1-1+width)%width][posY1]);
                H+=calcPairEnthalpy(ID1,map[posX1][(posY1-1+height)%height]);
    }
    #ifndef NDEBUG
    std::cerr<<"H "<<H<<std::endl;
    #endif
    
    return H;
}

double Lipidsystem::calcHostFreeEnerg()
{
    #ifndef NDEBUG
    std::cerr<<"Lipidsystem::calcHostFreeEnerg"<<std::endl;
    std::cerr<<"ID0 "<<ID0<<" posX0 "<<posX0<<" posY0 "<<posY0<<std::endl;
    #endif

    
    double G=0;

    G+=calcPairEnthalpy(ID0,map[(posX0-1+width)%width][posY0]);
    G+=calcPairEnthalpy(ID0,map[(posX0+1)%width][posY0]);
    G+=calcPairEnthalpy(ID0,map[posX0][(posY0+1)%height]);
    G+=calcPairEnthalpy(ID0,map[posX0][(posY0-1+height)%height]);

   
    #ifndef NDEBUG
    std::cerr<<"H "<<G<<std::endl;
    std::cerr<<"kB T S "<<-inputfile->kBT*lipidproperties->entropyFunction[lipids[ID0].getType()][lipids[ID0].getOrderPara()]<<std::endl;
    std::cerr<<"self E "<<lipidproperties->selfEnergieFunction[lipids[ID0].getType()][lipids[ID0].getOrderPara()]<<std::endl;
    #endif

    G+=lipidproperties->selfEnergieFunction[lipids[ID0].getType()][lipids[ID0].getOrderPara()]-inputfile->kBT*lipidproperties->entropyFunction[lipids[ID0].getType()][lipids[ID0].getOrderPara()];
    
    
    return G;
}




double Lipidsystem::calcPairEnthalpy(unsigned int ID_1,unsigned int ID_2)
{ 
    int type1=lipids[ID_1].getType();
    int type2=lipids[ID_2].getType();
    int order1=lipids[ID_1].getOrderPara();
    int order2=lipids[ID_2].getOrderPara();
    

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
    std::cerr<<"Lipidsystem::fluctuate"<<std::endl;
    #endif

    oldOrder=lipids[ID0].getOrderPara();
    int maxOrder=inputfile->types[lipids[ID0].getType()].maxOrder;
    int minOrder=inputfile->types[lipids[ID0].getType()].minOrder;
    int maxFluc=inputfile->types[lipids[ID0].getType()].maxFluc;
    
    lipids[ID0].setOrderPara(enhance::random_int((oldOrder-maxFluc+minOrder+std::abs(oldOrder-maxFluc-minOrder))/2,(oldOrder+maxFluc+maxOrder-std::abs(oldOrder+maxFluc-maxOrder))/2));
    
    #ifndef NDEBUG
    std::cerr<<"oldOrder: "<<oldOrder<<" maxOrder: "<<maxOrder<<" minOrder: "<<minOrder<<" maxFluc: "<<maxFluc<<std::endl;
    std::cerr<<"searching order between "<<(oldOrder-maxFluc+minOrder+std::abs(oldOrder-maxFluc-minOrder))/2<<" "<<(oldOrder+maxFluc+maxOrder-std::abs(oldOrder+maxFluc-maxOrder))/2<<" new order: "<<lipids[ID0].getOrderPara()<<std::endl;
    #endif
    
    
    
}

void Lipidsystem::fluctuateBack()
{
    #ifndef NDEBUG
    std::cerr<<"Lipidsystem::fluctuateBack"<<std::endl;
    #endif
    lipids[ID0].setOrderPara(oldOrder);
    
}







