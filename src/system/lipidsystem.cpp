#include "lipidsystem.h"

Lipidsystem::Lipidsystem()
{

}

void Lipidsystem::readParas( std::shared_ptr<LipidProperties> _lipidproperties,std::shared_ptr<InputFile> _inputfile)
{
    lipidproperties=_lipidproperties;
    inputfile=_inputfile;
    
    map =new  int*[inputfile->width];
    for(int i=0;i<inputfile->width;i++)
        map[i] =new  int[inputfile->height];
}


void Lipidsystem::setup()
{   
    #ifndef NDEBUG
    std::cerr<<"Lipidsystem::setup"<<std::endl;
    #endif
    
    
    int defaultOrderIndex=(inputfile->paras.at("defaultOrderPara")-inputfile->paras.at("minOrder"))/inputfile->paras.at("DeltaOrder");
    
    //setting up type distr
    int N=inputfile->width*inputfile->height;
    std::vector<int> types{};
    for(int type=0;type<inputfile->nType-1;type++)
        for(int j=0;j<inputfile->concentrations[type]*N and types.size()<=N ;j++)
            types.push_back(type);

    while (types.size()<N) types.push_back(inputfile->nType-1); 
    std::shuffle(types.begin(), types.end(), enhance::rand_engine); //randomize order of types


    for(int i=0;i<inputfile->width;i++)
    {
        for(int j=0;j<inputfile->height;j++)
        {   
            map[i][j]=i*inputfile->height+j;
            lipids.push_back(Lipid(i*inputfile->height+j,types[i*inputfile->height+j],defaultOrderIndex,i,j));
        }
    }
}


void Lipidsystem::setOrder(boost::multi_array<int, 2> data)
{
    #ifndef NDEBUG
    std::cout<<"Lipidsystem::setOrder"<<std::endl;
    #endif
    for(int i=0;i<inputfile->width;i++)
        for(int j=0;j<inputfile->height;j++)
            lipids[i*inputfile->height+j].setOrder(data[i][j]);

}

void Lipidsystem::setTypes(boost::multi_array<int, 2> data)
{
    #ifndef NDEBUG
    std::cout<<"Lipidsystem::setTypes"<<std::endl;
    #endif
    for(int i=0;i<inputfile->width;i++)
        for(int j=0;j<inputfile->height;j++)
            lipids[i*inputfile->height+j].setType(data[i][j]);

}


int Lipidsystem::getMeanOrder()
{
    #ifndef NDEBUG
    std::cerr<<"Lipidsystem::getMeanOrder"<<std::endl;
    #endif
    
    int mean=0;
    int count=0;
    for(int i=0;i<inputfile->width;i++)
    {
        for(int j=0;j<inputfile->height;j++)
        {   
            count++;
            mean += lipids[map[i][j]].getOrder();
        }
    }
    return mean/count;
}

std::vector<int> Lipidsystem::getOrderDistr()
{
    std::vector<int> distr((int)inputfile->paras.at("maxOrderIndex")+1,0);
    for(int i=0;i<inputfile->width;i++)
    for(int j=0;j<inputfile->height;j++)
        distr[lipids[map[i][j]].getOrder()]++;
    return distr;
    
}


const boost::multi_array<int,2> Lipidsystem::getOrderParas()
{
    boost::multi_array<int,2> data(boost::extents[inputfile->width][inputfile->height]);
    for(int i=0;i<inputfile->width;i++)    
    for(int j=0;j<inputfile->height;j++)
        data[i][j] = lipids[map[i][j]].getOrder();
    return data;
}

const boost::multi_array<int,2> Lipidsystem::getTypes()
{
    boost::multi_array<int,2> data(boost::extents[inputfile->width][inputfile->height]);
    for(int i=0;i<inputfile->width;i++)
    for(int j=0;j<inputfile->height;j++)  
        data[i][j] = lipids[map[i][j]].getType();
    return data;
}

const boost::multi_array<int,2> Lipidsystem::getIDs()
{
    boost::multi_array<int,2> data(boost::extents[inputfile->width][inputfile->height]);
    for(int i=0;i<inputfile->width;i++)
    for(int j=0;j<inputfile->height;j++)  
        data[i][j] = map[i][j];
    return data;
}




void Lipidsystem::printMap()
{
    for(int i=0;i<inputfile->width;i++)
    {
        for(int j=0;j<inputfile->height;j++)
            std::cout<<map[i][j]<<" ";
        std::cout<<std::endl;
    }
}

void Lipidsystem::swap(int ID_0, int ID_1)
{
    #ifndef NDEBUG
    std::cerr<<"Lipidsystem::swap"<<std::endl;
    #endif

    std::swap(lipids[ID_0].posX,lipids[ID_1].posX);
    std::swap(lipids[ID_0].posY,lipids[ID_1].posY);
    std::swap(map[lipids[ID_0].posX][lipids[ID_0].posY],map[lipids[ID_1].posX][lipids[ID_1].posY]);
}


void Lipidsystem::fluctuate(int ID)
{
    #ifndef NDEBUG
    std::cerr<<"Lipidsystem::fluctuate"<<std::endl;
    #endif

    oldOrder=lipids[ID].getOrder();
    int maxOrder=std::get<1>(inputfile->types[lipids[ID].getType()]);
    int minOrder=std::get<2>(inputfile->types[lipids[ID].getType()]);
    int maxFluc=std::get<3>(inputfile->types[lipids[ID].getType()]);
    
    lipids[ID].setOrder(enhance::random_int((oldOrder-maxFluc+minOrder+std::abs(oldOrder-maxFluc-minOrder))/2,(oldOrder+maxFluc+maxOrder-std::abs(oldOrder+maxFluc-maxOrder))/2));
    
    #ifndef NDEBUG
    std::cerr<<"oldOrder: "<<oldOrder<<" maxOrder: "<<maxOrder<<" minOrder: "<<minOrder<<" maxFluc: "<<maxFluc<<std::endl;
    std::cerr<<"searching order between "<<(oldOrder-maxFluc+minOrder+std::abs(oldOrder-maxFluc-minOrder))/2<<" "<<(oldOrder+maxFluc+maxOrder-std::abs(oldOrder+maxFluc-maxOrder))/2<<" new order: "<<lipids[ID].getOrder()<<std::endl;
    #endif
    
    
    
}

void Lipidsystem::fluctuateBack(int ID)
{
    #ifndef NDEBUG
    std::cerr<<"Lipidsystem::fluctuateBack"<<std::endl;
    #endif
    lipids[ID].setOrder(oldOrder);
    
}







