#ifndef TYPEPROPERTIES_H
#define TYPEPROPERTIES_H
#include <iostream>

class TypeProperties
{  
    public:
        TypeProperties(std::string,int,int,int);
        std::string typeName;
        int maxOrder=0;
        int minOrder=0;
        int maxFluc;
//         int getOrderIndex(double&);
//         double getOrderValue(int&);

};

#endif // TYPEPROPERTIES_H
