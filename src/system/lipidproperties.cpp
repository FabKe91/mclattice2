#include "lipidproperties.h"


LipidProperties::LipidProperties()
{
}


double NN_DPPC(double temp, double order)
{
   double Tm        = 323;
   double mag       = 10.;
   double sign      = -2*(Tm<temp)+1;
   double x_Tshift  = 0.65 + 0.4 / (1 + std::exp( mag * (-(temp - Tm)) ) );
   double y_Tc      = 0.61 - 0.7 * sign;
   double y_Tf      = 0.042 + 0.018 * sign;
   double y_Tmag    = 1.5 - 1.0 * sign;
   double y_Tshift  = sign * (y_Tmag / (1 + std::exp(-(y_Tf * std::abs(Tm - temp)))) ) ;
   double k         = 1.8 + 12.2 / (1 + std::exp( mag * (-(Tm - temp)) ) );

   return ( 4.12 + y_Tc + y_Tshift + sign * 0.5 / (1 + std::exp(-k * (order - x_Tshift))) );
}



// double NN_DUPC(double temp, double order)
// {
//    return 4.5943514+0.0370882*order-0.0037423*temp;
// }


void LipidProperties::readParas(std::shared_ptr<InputFile> _inputfile)
{
    #ifndef NDEBUG
    std::cout<<"LipidProperties::readParas"<<std::endl;
    #endif

    inputfile=_inputfile;

   
    //array construction
    neighbourFunction= new double*[inputfile->nType];
    entropyFunction= new double*[inputfile->nType];
    selfEnergieFunction= new double*[inputfile->nType];
    lipidCholEnergieFunction= new double**[inputfile->nType];
    enthalpyFunction= new double***[inputfile->nType];
    cholLipidNeigh= new double[5];
    cholCholEnergie= new double[5];
    
    
    
    for(int i=0;i<inputfile->nType;i++)
    {
        neighbourFunction[i]= new double[(int)inputfile->paras.at("maxOrderIndex")+1];
        entropyFunction[i]= new double[(int)inputfile->paras.at("maxOrderIndex")+1];
        selfEnergieFunction[i]= new double[(int)inputfile->paras.at("maxOrderIndex")+1];
        lipidCholEnergieFunction[i]= new double*[6];
        enthalpyFunction[i]= new double**[i+1];
        
        
        for(int j=0;j<6;j++)
            lipidCholEnergieFunction[i][j]= new double[(int)inputfile->paras.at("maxOrderIndex")+1];



        for(int j=0;j<=i;j++)
        {  
            enthalpyFunction[i][j]=new double* [7];
            for(int l=0;l<7;l++)
                enthalpyFunction[i][j][l]=new double [(int)inputfile->paras["maxOrderIndex"]+1];
        }
    }
    
    
    //set array values
    
    for(int i=0;i<5;i++)
    {
        cholLipidNeigh[i]=enhance::polynom(inputfile->cholLipidNeighPara[i],inputfile->paras.at("T"));
//         std::cout<<cholLipidNeigh[i]<<std::endl;
        cholCholEnergie[i]=enhance::polynom(inputfile->CholCholEnergiePara,i);
    }
    
    
    
//     for(int i=0;i<inputfile->nType;i++)
//     {
//         for(int j=0;j<5;j++)
//         {
//             neighbourFunction[i][j]=enhance::polynom(inputfile->LipidLipidNeighPara[i][j],inputfile->paras.at("T"));
//         }
//     }
//     
    
    for(int i=0;i<inputfile->nType;i++)
    {
            int k=0;

        for(double order=inputfile->paras.at("minOrder");order<inputfile->paras.at("maxOrder")+inputfile->paras.at("DeltaOrder");order+=inputfile->paras.at("DeltaOrder"))
        {   
            
            if (std::get<0>(inputfile->types[i])=="DPPC")   neighbourFunction[i][k]=NN_DPPC(inputfile->paras.at("T"),order);
//             else if (std::get<0>(inputfile->types[i])=="DUPC")   neighbourFunction[i][k]=NN_DUPC(inputfile->paras.at("T"),order);
            else throw std::invalid_argument("no NN funktion found for type: "+std::get<0>(inputfile->types[i]));

            entropyFunction[i][k]=enhance::polynom(inputfile->entropyPara[i],order);
            selfEnergieFunction[i][k]=enhance::polynom(inputfile->selfEnergiePara[i],order);
            
            for(int j=0;j<6;j++)
                lipidCholEnergieFunction[i][j][k]=enhance::polynom(inputfile->lipidCholEnergiePara[i][j],order);

            for(int j=0;j<=i;j++)
            {
                for(int l=0;l<7;l++)
                    enthalpyFunction[i][j][l][k]=enhance::polynom(inputfile->enthalpyPara[i][j][l],order);
            }
        
            k++;
        }
    }
}



// void LipidProperties::updateKBT() //not used currently
// {
//     for(int i=0;i<inputfile->nType;i++)
//     {
//             int k=0;
//         for(double order=inputfile->paras.at("minOrder");order<inputfile->paras.at("maxOrder")+inputfile->paras.at("DeltaOrder");order+=inputfile->paras.at("DeltaOrder"))
//         {   
//             if (inputfile->types[i].typeName=="DPPC")   neighbourFunction[i][k]=NN_DPPC(inputfile->paras.at("T"),order);
//             else if (inputfile->types[i].typeName=="DUPC")   neighbourFunction[i][k]=NN_DUPC(inputfile->paras.at("T"),order);
//             else throw std::invalid_argument("no NN funktion found for type: "+inputfile->types[i].typeName);
//             k++;
//         }
//     }
// 
//     
// }




