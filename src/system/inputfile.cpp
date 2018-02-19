#include "inputfile.h"

InputFile::InputFile(std::string filename)
{
    std::ifstream infile(filename);
    std::string line;
    double val;  //buffer
    std::string name;
    std::string name2;
    double DeltaOrder;
    bool firstCall =true; //see name=="Enthalpy"

    while (std::getline(infile, line))
    {   
        
        
        //erase comments
        auto pos=line.find("#");
        if (pos==0 or line.size()==0) continue;
        if (pos==std::string::npos);
        else
        {
            line.erase(pos,std::string::npos);
        }
//         std::cout<<line<<std::endl;

        std::istringstream iss(line);
        if(!(iss>>name)) throw std::invalid_argument( "cant read line: "+line );
        if(name=="outs")
        {
            while(iss>>name)
            {
                outs.push_back(name);
            }
        }
        else if(name=="Type") //type definition 
        {
            std::string typeName;
            double maxOrder;
            double minOrder;
            double maxFluc;

            iss>>typeName>>name>>minOrder>>name>>maxOrder>>name>>maxFluc; //discarding parameter names, only order relevant
            types.push_back(TypeProperties(typeName,(maxOrder-paras["minOrder"])/paras["DeltaOrder"],(minOrder-paras["minOrder"])/paras["DeltaOrder"],maxFluc/paras["DeltaOrder"]));
            typeMap[typeName]=nType;
            
            std::vector<double> vec;
            neighbourPara.push_back(vec);
            entropyPara.push_back(vec);
            nType++;
        }
        else if(name=="Neighbour")
        {
            iss>>name;
            while(iss>>val)
            {
            neighbourPara[typeMap[name]].push_back(val);
            }
                
        }        
        else if(name=="Entropy")
        {
            iss>>name;
            while(iss>>val)
            {
            entropyPara[typeMap[name]].push_back(val);
            }
                
        }
        else if(name=="Enthalpy")
        {
            if (firstCall) //if first called, create vec<vec<>>... for defined types to push_back fit parameters
            {
                std::vector<double> vec;
                std::vector<std::vector<double>> vec2;
                for(int i=0;i<nType;i++)
                {
                    vec2.push_back(vec);
                    enthalpyPara.push_back(vec2);
                }
                firstCall=false;
            }
            iss>>name>>name2;
            while(iss>>val)
            {
                if (typeMap[name]>=typeMap[name2])
                {
                    enthalpyPara[typeMap[name]][typeMap[name2]].push_back(val);

                }                
                else
                {
                    enthalpyPara[typeMap[name2]][typeMap[name]].push_back(val);
                }
            }
        }
        else
        {
            if(!(iss>>val)) throw std::invalid_argument( "can't read line: "+line );
            paras[name]=val;    
//             std::cout<<name<<paras[name]<<std::endl;
        }
    }
    
    
    paras["maxOrderIndex"]=(int)((paras["maxOrder"]-paras["minOrder"])/paras["DeltaOrder"]); //set max Order [-0.5,1]->[0,150]
    

    //check for parameters
    for(int i=0;i<nType;i++)
    {
        if(neighbourPara[i].size()==0) throw std::invalid_argument( "missing neighbourPara" );        
        if(entropyPara[i].size()==0) std::cout<<"no entropy given, can't run jet"<<std::endl;
        for(int j=0;j<=i;j++)
        {
            if(enthalpyPara[i][j].size()==0) throw std::invalid_argument( "missing enthalpyPara" );
        }
        
    }
    
}
