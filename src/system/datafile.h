#ifndef DATAFILE_H
#define DATAFILE_H
#include "hdf5.h"
#include <string>
#include <vector>

#include <iostream>
#include "inputfile.h"
#include <boost/multi_array.hpp>
#include "lipidsystem.h"
#include <iostream>
#include <functional>
#include <vector>
#include <map>
#include <stdio.h>

class DataFile
{

private:    
    const char* filename="out.h5";
    
    std::vector<std::string> outs;
    std::map<std::string,std::function<const boost::multi_array<int,2>()>> getterMap;
    std::vector<boost::multi_array<int,3>> buffer;
    int bufferLen=0;
    long int bufferSize=0;
    long int maxBufferSize=1024*1024*10;
    
    
    unsigned int   NX = 0;          
    unsigned int   NY = 0;
    unsigned int images = 0;

    hsize_t     dimsf[3];
    hsize_t     offset[3]={0,0,0};
    hsize_t      size[3];
    hsize_t      chunk_dims[3];
    hsize_t      maxdimsf[3] = {H5S_UNLIMITED, H5S_UNLIMITED,H5S_UNLIMITED};

    hid_t spaceDummy;
    std::shared_ptr<InputFile> inputfile;

    
public:


    ~DataFile();
    DataFile(Lipidsystem&,std::shared_ptr<InputFile>);
    void createFile();
    void writeStep();

private:
    template<typename INTorFloat>
    void extendDataset(std::string ,const boost::multi_array<INTorFloat,3>&,hid_t&);
    void createDataset(std::string, hid_t&);
    void createAttribute(std::string, double, hid_t&);
    void flush();
    
};

#endif // DATAFILE_H
