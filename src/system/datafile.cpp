#include "datafile.h"




DataFile::~DataFile()
{
    flush();
}


DataFile::DataFile(Lipidsystem& lipidsystem,std::shared_ptr<InputFile> _inputfile)
{
    inputfile=_inputfile;
    
    getterMap["orderPara"]=std::bind(&Lipidsystem::getOrderParas,&lipidsystem);
    getterMap["Type"]=std::bind(&Lipidsystem::getTypes,&lipidsystem);
    getterMap["IDs"]=std::bind(&Lipidsystem::getIDs,&lipidsystem);


    NX=inputfile->width;
    NY=inputfile->height;
    maxBufferSize=1024*1024*inputfile->paras.at("maxBufferSize");

}


void DataFile::createFile()
{
    hid_t file;
    
    file= H5Fcreate( filename, H5F_ACC_TRUNC , H5P_DEFAULT, H5P_DEFAULT);
    
    dimsf[0] = 1;
    dimsf[1] = NX;
    dimsf[2] = NY;

    size[0] = 0;
    size[1] = NX;
    size[2] = NY;

    chunk_dims[0] = 1;
    chunk_dims[1] = NX;
    chunk_dims[2] = NY;


    for(auto& o: inputfile->outs)
    {
        createDataset(o,file);
    }
    
//     for (auto const& mapitem: inputfile->paras) createAttribute(mapitem.first,mapitem.second,file);
    

}

void DataFile::createDataset(std::string datasetName, hid_t& file)
{
    hid_t   dataspace,dcpl,dataset ;
//     herr_t          status;

    dataspace= H5Screate_simple( 3, size, NULL);  

      
    dcpl = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk (dcpl, 3, chunk_dims);

   
    
    dataset = H5Dcreate (file, datasetName.c_str(), H5T_NATIVE_INT, dataspace, H5P_DEFAULT, dcpl,
                H5P_DEFAULT); 
    
    buffer.push_back(boost::multi_array<int,3>(boost::extents[0][NX][NY]));
}

// void DataFile::createAttribute(std::string attrName, double val, hid_t& file)
// {
//     double* data;
//     data=&val;
//     const hsize_t dims=1;
//     DataSpace* dspace = new DataSpace(1,&dims);
//     Attribute attr = file.createAttribute(attrName, PredType::NATIVE_DOUBLE, *dspace);
//     delete dspace;
//     
//     attr.write(PredType::NATIVE_DOUBLE,data);
//     
// }

void DataFile::writeStep()
{
    #ifndef NDEBUG
    std::cout<<"DataFile::writeStep"<<std::endl;
    #endif

    bufferLen++;
    for(int i=0; i<inputfile->outs.size();i++)
    {
        buffer[i].resize(boost::extents[bufferLen][NX][NY]);
        buffer[i][bufferLen-1]=getterMap.at(inputfile->outs[i])();

    }
    bufferSize+=NX*NY*sizeof(int)*inputfile->outs.size();
    if(bufferSize>=maxBufferSize) flush();

}

void DataFile::flush()
{
    hid_t file;
    
    file= H5Fcreate( filename, H5F_ACC_RDWR , H5P_DEFAULT, H5P_DEFAULT);
    
    images+=bufferLen;
    dimsf[0]=bufferLen;
    offset[0]=images-bufferLen;
    size[0]=images;
    
    spaceDummy=H5Screate_simple(3, dimsf ,NULL); 

    for(int i=0; i<inputfile->outs.size();i++)
    {
        extendDataset(inputfile->outs[i],buffer[i],file);
        buffer[i].resize(boost::extents[0][NX][NY]);
    }
    
    std::cout<<"flushing! size: "<<bufferSize/1024.0/1024<<"MB images: "<<images<<std::endl;

    bufferLen=0;
    bufferSize=0;
    
}


template<typename INTorFloat>
void DataFile::extendDataset(std::string datasetName, const boost::multi_array<INTorFloat,3>& data_array, hid_t& file)
{  
    hid_t dataset,filespace;
    dataset=H5Dopen(file, datasetName.c_str(),H5P_DEFAULT);
    
    
    H5Dset_extent(dataset,size);
    filespace = H5Dget_space(dataset);
    
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset,NULL,dimsf, NULL );
    H5Dwrite(dataset,  H5T_NATIVE_INT, spaceDummy, filespace, H5P_DEFAULT,data_array.data() );

}


