#include "inputfile.h"

InputFile::InputFile(std::string filename)
{
    #ifndef NDEBUG
    std::cout<<"InputFile::InputFile"<<std::endl;
    std::cout<<"reading lines..."<<std::endl;
    #endif
    
    std::ifstream infile(filename);
    std::string line;
    double val;  //buffer for values
    double val2;  //buffer for values
    std::string name; //buffer for strings
    bool firstCall =true; //see name=="Enthalpy"

    while (std::getline(infile, line))
    {   
        #ifndef NDEBUG
        std::cout<<line<<std::endl;
        #endif      

        //remove comments
        auto pos=line.find("#");
        if (pos==0 or line.size()==0) continue;
        if (pos==std::string::npos);
        else
        {
            line.erase(pos,std::string::npos);
        }

        std::istringstream iss(line);
        if(!(iss>>name)) throw std::invalid_argument( "cant read line: "+line );
        if(name=="outs")    while(iss>>name)    outs.push_back(name);
        else if(name=="Type") //type definition 
        {
            if (!firstCall) throw std::invalid_argument( "type def must happen bevor Enthalpy "+line );
            std::string typeName;
            iss>>typeName;
            typeMap[typeName]=nType;
            std::map<std::string,double> typeParas;

            while(iss>>name>>val)   typeParas[name]=val;
            
            types.push_back(TypeProperties(typeName,(typeParas.at("maxOrder")-paras.at("minOrder"))/paras.at("DeltaOrder"),(typeParas.at("minOrder")-paras.at("minOrder"))/paras.at("DeltaOrder"),typeParas.at("maxFluc")/paras.at("DeltaOrder")));
            
            concentrations.push_back(typeParas.at("konz"));
            
            
            std::vector<double> vec;
            std::vector<std::vector<double>> vec2;
            
            LipidLipidNeighPara.push_back(vec2);
            entropyPara.push_back(vec);
            selfEnergiePara.push_back(vec);
            lipidCholEnergiePara.push_back(vec2);
            nType++;
        }
        else if(name=="LipidLipidNeighPara")
        {
            iss>>name;
            iss>>val;
            
            std::vector<double> vec;
            while(iss>>val2)
            {
                vec.push_back(val2);
            }
            LipidLipidNeighPara[typeMap.at(name)].push_back(vec);
        }        
        else if(name=="Entropy")
        {
            iss>>name;
            while(iss>>val)
            {
            entropyPara[typeMap.at(name)].push_back(val);
            }
                
        }
        else if(name=="SelfEnergie")
        {
            iss>>name;
            while(iss>>val)
            {
            selfEnergiePara[typeMap.at(name)].push_back(val);
            }
                
        }
        else if(name=="CholCholEnergie")
        {
            while(iss>>val)
            {
                CholCholEnergiePara.push_back(val);
            }
        }        
        else if(name=="LipidCholEnergie")
        {
            iss>>name;
            iss>>val;
            
            std::vector<double> vec;
            while(iss>>val2)
            {
                vec.push_back(val2);
            }
            lipidCholEnergiePara[typeMap.at(name)].push_back(vec);

        }
        else if(name=="CholLipidNeigh")
        {
            std::vector<double> vec;
            cholLipidNeighPara.push_back(vec);
            iss>>val;
            while(iss>>val2)
            {
            cholLipidNeighPara[val].push_back(val2);
            }
                
        }
        else if(name=="Enthalpy")
        {
            if (firstCall) //if first called, create vec<vec<>>... for defined types to push_back fit parameters
            {
                std::vector<std::vector<double>> vec;
                std::vector<std::vector<std::vector<double>>> vec2;
                for(int i=0;i<nType;i++)
                {
                    vec2.push_back(vec);
                    enthalpyPara.push_back(vec2);
                }
                firstCall=false;
            }
            std::string name2;
            iss>>name>>name2;
            iss>>val;
            std::vector<double> vec3;
            while(iss>>val)
            {
                vec3.push_back(val);
            }
            if (typeMap.at(name)>=typeMap.at(name2))
            {
                enthalpyPara[typeMap.at(name)][typeMap.at(name2)].push_back(vec3);

            }                
            else
            {
                enthalpyPara[typeMap.at(name2)][typeMap.at(name)].push_back(vec3);
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
    
    T=paras.at("T");
    kBT=paras.at("kB")*T;
    width=paras.at("width");
    height=paras.at("height");
    
    paras["numberCholesterin"]=paras.at("height")*paras.at("width")/(1/paras.at("cholesterinConc")-1);
    std::cout<<paras["numberCholesterin"]<<std::endl;
    
    
    
    
    

    //check for parameters
    for(int i=0;i<nType;i++)
    {
//         if(neighbourPara[i].size()==0) throw std::invalid_argument( "missing neighbourPara" );        
        if(entropyPara[i].size()==0) throw std::invalid_argument( "no entropy given, can't run jet" );
        if(selfEnergiePara[i].size()==0) throw std::invalid_argument("no selfEnergie given, can't run jet");
        for(int j=0;j<=i;j++)
        {
            if(enthalpyPara[i][j].size()==0) throw std::invalid_argument( "missing enthalpyPara" );
        }
        
    }
    #ifndef NDEBUG
    std::cout<<"InputFile::InputFile   finished reading"<<std::endl;
    #endif
    
}
