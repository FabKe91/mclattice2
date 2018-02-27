#include "lipidproperties.h"


LipidProperties::LipidProperties()
{
}


double NN_DPPC(double temp, double order)
{
   double Tm        = 323;
   double mag       = 10.;
   double sign      = -2*Tm<temp+1;
   double x_Tshift  = 0.65 + 0.4 / (1 + std::exp( mag * (-(temp - Tm)) ) );
   double y_Tc      = 0.61 - 0.7 * sign;
   double y_Tf      = 0.042 + 0.018 * sign;
   double y_Tmag    = 1.5 - 1.0 * sign;
   double y_Tshift  = sign * (y_Tmag / (1 + std::exp(-(y_Tf * std::abs(Tm - temp)))) ) ;
   double k         = 1.8 + 12.2 / (1 + std::exp( mag * (-(Tm - temp)) ) );

   return ( 4.12 + y_Tc + y_Tshift + sign * 0.5 / (1 + std::exp(-k * (order - x_Tshift))) );
}




double NN_DUPC(double temp, double order)
{
   return 4.5943514+0.0370882*order-0.0037423*temp;
}

void LipidProperties::readParas()
{
    #ifndef NDEBUG
    std::cout<<"LipidProperties::readParas"<<std::endl;
    #endif


        std::cout<<"maxOrderIndex"<<std::endl;
        std::cout<<InputFile::paras.at("maxOrderIndex")<<std::endl;

    
    //array construction
    neighbourFunction= new double*[InputFile::nType];
    entropyFunction= new double*[InputFile::nType];
    selfEnergieFunction= new double*[InputFile::nType];
    enthalpyFunction= new double**[InputFile::nType];
    
    for(int i=0;i<InputFile::nType;i++)
    {
        neighbourFunction[i]= new double[(int)InputFile::paras.at("maxOrderIndex")+1];
        entropyFunction[i]= new double[(int)InputFile::paras.at("maxOrderIndex")+1];
        selfEnergieFunction[i]= new double[(int)InputFile::paras.at("maxOrderIndex")+1];
        enthalpyFunction[i]= new double*[i+1];
        for(int j=0;j<=i;j++)
        {
            enthalpyFunction[i][j]=new double[(int)InputFile::paras["maxOrderIndex"]+1];
        }

    }
    
    
    //set array values
    for(int i=0;i<InputFile::nType;i++)
    {
            int k=0;

        for(double order=InputFile::paras.at("minOrder");order<InputFile::paras.at("maxOrder")+InputFile::paras.at("DeltaOrder");order+=InputFile::paras.at("DeltaOrder"))
        {   
            if (InputFile::types[i].typeName=="DPPC")   neighbourFunction[i][k]=NN_DPPC(InputFile::paras.at("T"),order);
            else if (InputFile::types[i].typeName=="DUPC")   neighbourFunction[i][k]=NN_DUPC(InputFile::paras.at("T"),order);
            else throw std::invalid_argument("no NN funktion found");
            
            
//             if (InputFile::types[i].typeName=="DPPC") neighbourFunction[i][k]=enhance::sigmoid(InputFile::neighbourPara[i],order);
//             else if (InputFile::types[i].typeName=="DUPC") neighbourFunction[i][k]=enhance::polynom(InputFile::neighbourPara[i],order);
            entropyFunction[i][k]=enhance::polynom(InputFile::entropyPara[i],order);
            selfEnergieFunction[i][k]=enhance::polynom(InputFile::selfEnergiePara[i],order);
//             std::cout<<selfEnergieFunction[i][k]<<", ";

            
            for(int j=0;j<=i;j++)
            {
                enthalpyFunction[i][j][k]=enhance::polynom(InputFile::enthalpyPara[i][j],order);
            }
        
            k++;
        }
    }
}



void LipidProperties::updateKBT()
{
    for(int i=0;i<InputFile::nType;i++)
    {
            int k=0;
        for(double order=InputFile::paras.at("minOrder");order<InputFile::paras.at("maxOrder")+InputFile::paras.at("DeltaOrder");order+=InputFile::paras.at("DeltaOrder"))
        {   
            if (InputFile::types[i].typeName=="DPPC")   neighbourFunction[i][k]=NN_DPPC(InputFile::paras.at("T"),order);
            else if (InputFile::types[i].typeName=="DUPC")   neighbourFunction[i][k]=NN_DUPC(InputFile::paras.at("T"),order);
            else throw std::invalid_argument("no NN funktion found");
            k++;
        }
    }

    
}




