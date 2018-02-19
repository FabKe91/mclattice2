#include "typeproperties.h"

TypeProperties::TypeProperties(std::string _typeName, int _maxOrder, int _minOrder, int _maxFluc)
{
    typeName=_typeName;
    minOrder=_minOrder;
    maxOrder=_maxOrder;
    maxFluc=_maxFluc;
//     std::cout<<minOrder<<maxOrder<<std::endl;
    
}

// int TypeProperties::getOrderIndex(double& orderValue)
// {
//     return (orderValue-minOrderRaw)/DeltaOrder;
// }

// double TypeProperties::getOrderValue(int& orderIndex)
// {
//     return orderIndex*DeltaOrder+minOrderRaw;
// }
