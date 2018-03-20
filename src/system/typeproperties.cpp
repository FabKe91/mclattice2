#include "typeproperties.h"

TypeProperties::TypeProperties(std::string _typeName, int _maxOrder, int _minOrder, int _maxFluc)
{
    typeName=_typeName;
    minOrder=_minOrder;
    maxOrder=_maxOrder;
    maxFluc=_maxFluc;    
}

