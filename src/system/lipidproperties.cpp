#include "lipidproperties.h"


LipidProperties::LipidProperties()
{
}


double DavitNN(double temp, double order)
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

void LipidProperties::readParas(std::shared_ptr<InputFile> _inputfile)
{
    #ifndef NDEBUG
    std::cout<<"LipidProperties::readParas"<<std::endl;
    #endif

    inputfile=_inputfile;

    
    
    //array construction
    neighbourFunction= new double*[inputfile->nType];
    entropyFunction= new double*[inputfile->nType];
    enthalpyFunction= new double**[inputfile->nType];
    
    for(int i=0;i<inputfile->nType;i++)
    {
        neighbourFunction[i]= new double[(int)inputfile->paras["maxOrderIndex"]+1];
        entropyFunction[i]= new double[(int)inputfile->paras["maxOrderIndex"]+1];
        enthalpyFunction[i]= new double*[i+1];
        for(int j=0;j<=i;j++)
        {
            enthalpyFunction[i][j]=new double[(int)inputfile->paras["maxOrderIndex"]+1];
        }

    }
    
    
    //set array values
    for(int i=0;i<inputfile->nType;i++)
    {
            int k=0;

        for(double order=inputfile->paras["minOrder"];order<inputfile->paras["maxOrder"]+inputfile->paras["DeltaOrder"];order+=inputfile->paras["DeltaOrder"])
        {   
            
            neighbourFunction[i][k]=DavitNN(inputfile->paras["T"],order);
            entropyFunction[i][k]=polynom(inputfile->entropyPara[i],order);

            
            for(int j=0;j<=i;j++)
            {
                enthalpyFunction[i][j][k]=polynom(inputfile->enthalpyPara[i][j],order);
            }
        
            k++;
        }
    }
}





double LipidProperties::polynom(std::vector<double>& coeff, double x)
{
    double y=0;
    for(int n=0; n<coeff.size();n++)
    {
        y+=std::pow(x,n)*coeff[n];
    }
    return y;
}

double LipidProperties::sigmoid(std::vector<double>& coeff, double x)
{
    if(coeff.size()!=4) throw std::invalid_argument("need 4 parameters for sigmoid");
    return coeff[0]  / (1 + std::exp(-coeff[1]  * (x - coeff[2] ))) + coeff[3];  
}
