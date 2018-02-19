#ifndef LIPID_H
#define LIPID_H


class Lipid
{

private:
    int  ID {0};
    int   type {0};
    int   orderPara {0};

public:
    Lipid(int,int,int);


    int  getID()   const { return ID; }
    int  getType() const { return type; }
    int  getOrderPara() const { return orderPara; }
    void  setOrderPara(int newOrder) { orderPara=newOrder; }


};

#endif // LIPID_H
