#include "cholesterinsystem.h"


CholesterinSystem::CholesterinSystem()
{
}

void CholesterinSystem::setup(std::shared_ptr<InputFile> _inputfile)
{
    inputfile=_inputfile;
    map =new  int*[inputfile->width];
    for(int i=0;i<inputfile->width;i++)
        map[i] =new  int[inputfile->height];
    
    
    
    
    int N=inputfile->width*inputfile->height;
    std::vector<bool> setChol{};
    
    
    for(int i=0;i<N;i++)
    {
        if (i<inputfile->paras.at("numberCholesterin")) setChol.push_back(true);
        else setChol.push_back(false);
    }
    std::shuffle(setChol.begin(), setChol.end(), enhance::rand_engine); //randomize order of chols


    for(int i=0;i<inputfile->width;i++)
    {
        for(int j=0;j<inputfile->height;j++)
        {   
            map[i][j]=inputfile->height*i+j;
            chols.push_back(Cholesterin(inputfile->height*i+j,i,j,setChol[i*inputfile->height+j]));
        }
    }
}

void CholesterinSystem::swap(int ID0, int ID1)
{
    #ifndef NDEBUG
    std::cout<<"CholesterinSystem::swap"<<std::endl;
    #endif
    
//     std::cout<<"ID0 "<<ID0<<" ID1 "<<ID1<<std::endl;
//     std::cout<<"chols[ID0].posX "<<chols[ID0].posX<<" chols[ID1].posX "<<chols[ID1].posX<<std::endl;
//     std::cout<<"chols[ID0].posY "<<chols[ID0].posY<<" chols[ID1].posY "<<chols[ID1].posY<<std::endl;
    
    
//     printMap();
    
//     std::cout<<"swapped"<<std::endl;
    std::swap(chols[ID0].posX,chols[ID1].posX);
    std::swap(chols[ID0].posY,chols[ID1].posY);
    std::swap(map[chols[ID0].posX][chols[ID0].posY],map[chols[ID1].posX][chols[ID1].posY]);
    
//     printMap();
}

const boost::multi_array<int,2> CholesterinSystem::getIDs()
{
    boost::multi_array<int,2> data(boost::extents[inputfile->width][inputfile->height]);
    for(int i=0;i<inputfile->width;i++)
        for(int j=0;j<inputfile->height;j++)  
            data[i][j] = map[i][j];
    return data;
}

const boost::multi_array<int,2> CholesterinSystem::getOcc()
{
    boost::multi_array<int,2> data(boost::extents[inputfile->width][inputfile->height]);
    for(int i=0;i<inputfile->width;i++)
        for(int j=0;j<inputfile->height;j++)  
            data[i][j] = chols[map[i][j]].occupied;
    return data;
}

void CholesterinSystem::printMap()
{
    for(int i=0;i<inputfile->width;i++)
    {
        for(int j=0;j<inputfile->height;j++)
            std::cout<<map[i][j]<<" ";
        std::cout<<std::endl;
    }
}


void CholesterinSystem::setOccupation(boost::multi_array<int, 2> data)
{
    #ifndef NDEBUG
    std::cout<<"CholesterinSystem::setOccupation"<<std::endl;
    #endif
    for(int i=0;i<inputfile->width;i++)
        for(int j=0;j<inputfile->height;j++)
            chols[map[i][j]].occupied=data[i][j];
}






