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

    
    dimsf[0] = 1;
    dimsf[1] = NX;
    dimsf[2] = NY;

    size[0] = 0;
    size[1] = NX;
    size[2] = NY;

    chunk_dims[0] = 1;
    chunk_dims[1] = NX;
    chunk_dims[2] = NY;


}


void DataFile::createFile()
{
    hid_t file;
    
    file= H5Fcreate( filename, H5F_ACC_TRUNC , H5P_DEFAULT, H5P_DEFAULT);


    for(auto& o: inputfile->outs)
    {
        createDataset(o,file);
    }
    

    for (auto const& mapitem: inputfile->paras) createAttribute(mapitem.first,mapitem.second,file);
    
    H5Fclose (file);
}

void DataFile::createDataset(std::string datasetName, hid_t& file)
{
    hid_t   dataspace,dcpl,dataset ;
//     herr_t          status;

    dataspace= H5Screate_simple( 3, size, maxdimsf);  

      
    dcpl = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_deflate (dcpl, 9); // 9 = best compression, lowest speed
    H5Pset_chunk (dcpl, 3, chunk_dims);

   
    
    dataset = H5Dcreate (file, datasetName.c_str(), H5T_NATIVE_INT, dataspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    
    
    H5Pclose (dcpl);
    H5Dclose (dataset);
    H5Sclose (dataspace);

    
    
    buffer.push_back(boost::multi_array<int,3>(boost::extents[0][NX][NY]));
}

void DataFile::createAttribute(std::string attrName, double val, hid_t& file)
{
    
    hid_t   space,attr;   


    space = H5Screate(H5S_SCALAR);
    
    attr = H5Acreate (file, attrName.c_str(), H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite (attr, H5T_NATIVE_DOUBLE, &val);
    

    H5Aclose (attr);
    H5Sclose (space);
}

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
    #ifndef NDEBUG
    std::cout<<"DataFile::flush"<<std::endl;
    #endif
    
    
    hid_t file;
    
    file= H5Fopen( filename, H5F_ACC_RDWR ,  H5P_DEFAULT);
    
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

    H5Fclose (file);
}


template<typename INTorFloat>
void DataFile::extendDataset(std::string datasetName, const boost::multi_array<INTorFloat,3>& data_array, hid_t& file)
{  
    #ifndef NDEBUG
    std::cout<<"DataFile::extendDataset"<<std::endl;
    #endif


    hid_t dataset,filespace;
    dataset=H5Dopen(file, datasetName.c_str(),H5P_DEFAULT);
    

    H5Dset_extent(dataset,size);

    filespace = H5Dget_space(dataset);
    spaceDummy=H5Screate_simple(3, dimsf ,NULL);


    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset,NULL,dimsf, NULL );
    H5Dwrite(dataset,  H5T_NATIVE_INT, spaceDummy, filespace, H5P_DEFAULT,data_array.data() );
    

    H5Dclose (dataset);
    H5Sclose (filespace);
}



boost::multi_array<int,2> DataFile::getLastStep(std::string datasetName)
{
    boost::multi_array<int,2> data(boost::extents[NX][NY]);
    hid_t file,dataset,filespace;
    
    
    file = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset = H5Dopen (file, datasetName.c_str(), H5P_DEFAULT);

    filespace = H5Dget_space(dataset);
    spaceDummy=H5Screate_simple(3, dimsf ,NULL);
    
    offset[0]=images-1;

    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset,NULL,dimsf, NULL );


    H5Dread (dataset, H5T_NATIVE_INT, spaceDummy, filespace, H5P_DEFAULT, data.data());
   

    H5Dclose (dataset);
    H5Sclose (filespace);
    H5Fclose (file);
    return data;
}


void DataFile::readFile()
{
    hid_t file,dataset,filespace;
    
    
    file = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset = H5Dopen (file, "orderPara", H5P_DEFAULT);

    filespace = H5Dget_space(dataset);


    hsize_t dims[3];
    H5Sget_simple_extent_dims(filespace, dims, NULL);
    images=dims[0];

    H5Dclose (dataset);
    H5Sclose (filespace);
    H5Fclose (file);
    
    
    for(auto& o: inputfile->outs)
    {
        buffer.push_back(boost::multi_array<int,3>(boost::extents[0][NX][NY]));
    }
}

