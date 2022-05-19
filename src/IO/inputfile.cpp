#include "inputfile.h"

InputFile::InputFile()
{
    #ifndef NDEBUG
    std::cout<<"InputFile::InputFile"<<std::endl;
    std::cout<<"reading in.txt..."<<std::endl;
    #endif
    
    std::ifstream infile("in.txt");
    std::string line;
    double val;  //buffer for values
    double val2;  //buffer for values
    std::string string; //buffer for strings
    std::string string2; //buffer for strings


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
        if(!(iss>>string)) throw std::invalid_argument( "in.txt: cant read line: "+line );
        if(string=="outs")    while(iss>>string)    outs.push_back(string);
        else if(string=="Type") //type definition 
        {
            std::string typeName;
            iss>>typeName;
            typeMap[typeName]=nType;
            std::map<std::string,double> typeParas;

            while(iss>>string>>val)   typeParas[string]=val; //reading type specific parameters
            
            //tuple values: 0->typeName 1->maxOrder 2->minOrder 3->maxFluc
            std::tuple<std::string,int,int,int> typeProperties=std::make_tuple(typeName,(typeParas.at("maxOrder")-paras.at("minOrder"))/paras.at("DeltaOrder"),(typeParas.at("minOrder")-paras.at("minOrder"))/paras.at("DeltaOrder"),typeParas.at("maxFluc")/paras.at("DeltaOrder"));
            types.push_back(typeProperties);
            
            concentrations.push_back(typeParas.at("konz"));
            
            
            std::vector<double> vec;
            std::vector<std::vector<double>> vec2;
            
            LipidLipidNeighPara.push_back(vec2);
            cholLipidNeighPara.push_back(vec2);
            entropyPara.push_back(vec);
            selfEnergiePara.push_back(vec);
            lipidCholEnergiePara.push_back(vec2);
            nType++;
        }
        else
        {
            if(!(iss>>val)) throw std::invalid_argument( "in.txt: can't read line: "+line );
            paras[string]=val;    
        }
    }
    paras["maxOrderIndex"]=(int)((paras["maxOrder"]-paras["minOrder"])/paras["DeltaOrder"]); //set max Order [-0.5,1]->[0,150]
    
    T=paras.at("T");
    kBT=paras.at("kB")*T;
    width=paras.at("width");
    height=paras.at("height");
    
    paras["numberCholesterin"]=paras.at("height")*paras.at("width")/(1/paras.at("cholesterinConc")-1);
    
    

    
    #ifndef NDEBUG
    std::cout<<"reading fitParameters.prm..."<<std::endl;
    #endif
    
    std::ifstream parameterFile("fitParameters.prm");

    while (std::getline(parameterFile, line))
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
        if(!(iss>>string)) throw std::invalid_argument( "fitParameters.prm: cant read line: "+line );

        else if(string=="LipidLipidNeighPara")
        {
            iss>>string;
            iss>>val;
            
            std::vector<double> vec;
            while(iss>>val2)
            {
                vec.push_back(val2);
            }
            LipidLipidNeighPara[typeMap.at(string)].push_back(vec);
        }        
        else if(string=="Entropy")
        {
            iss>>string;
            while(iss>>val)
            {
            entropyPara[typeMap.at(string)].push_back(val);
            }
                
        }
        else if(string=="SelfEnergie")
        {
            iss>>string;
            while(iss>>val)
            {
            selfEnergiePara[typeMap.at(string)].push_back(val);
            }
                
        }
        else if(string=="CholCholEnergie")
        {
            while(iss>>val)
            {
                CholCholEnergiePara.push_back(val);
            }
        }        
        else if(string=="LipidCholEnergie")
        {
            iss>>string;
            iss>>val;
            
            std::vector<double> vec;
            while(iss>>val2)
            {
                vec.push_back(val2);
            }
            lipidCholEnergiePara[typeMap.at(string)].push_back(vec);

        }
        else if(string=="CholLipidNeigh")

        {
            iss>>string;
            iss>>val;
            
            std::vector<double> vec;
            while(iss>>val2)
            {
                vec.push_back(val2);
            }
            cholLipidNeighPara[typeMap.at(string)].push_back(vec);
        }        
        else if(string=="Enthalpy")
        {
            std::vector<std::vector<double>> vec;
            std::vector<std::vector<std::vector<double>>> vec2;
            for(int i=0;i<nType;i++)
            {
                vec2.push_back(vec);
                enthalpyPara.push_back(vec2);
            }

            iss>>string>>string2;
            iss>>val;
            std::vector<double> vec3;
            while(iss>>val)
            {
                vec3.push_back(val);
            }
            if (typeMap.at(string)>=typeMap.at(string2))
            {
                enthalpyPara[typeMap.at(string)][typeMap.at(string2)].push_back(vec3);

            }                
            else
            {
                enthalpyPara[typeMap.at(string2)][typeMap.at(string)].push_back(vec3);
            }
        }
        else
        {
            throw std::invalid_argument( "fitParameters.prm: can't read line: "+line );
        }
    }
    
    

    
    
    

    //check for parameters
    for(int i=0;i<nType;i++)
    {
//         if(neighbourPara[i].size()==0) throw std::invalid_argument( "missing neighbourPara" );        
        if(entropyPara[i].size()==0) throw std::invalid_argument( "Setup failed. No entropy given." );
        if(selfEnergiePara[i].size()==0) throw std::invalid_argument("Setup failed. No selfEnergie given.");
        for(int j=0;j<=i;j++)
        {
            if(enthalpyPara[i][j].size()==0) throw std::invalid_argument( "Setup failed. Missing enthalpyPara" );
        }
        
    }
    #ifndef NDEBUG
    std::cout<<"InputFile::InputFile   finished reading"<<std::endl;
    #endif
    
}
