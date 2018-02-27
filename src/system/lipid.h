#ifndef LIPID_H
#define LIPID_H


class Lipid
{

private:
    int ID;
    int type;
    int orderPara;


public:
    Lipid(int,int,int,int,int);
    
    int posX;
    int posY;

//     int getX() const { return posX; }
//     int getY() const { return posY; }
//     void setX(int X) { posX=X; }
//     void setY(int Y) { posY=Y; }
    int  getID()   const { return ID; }
    int  getType() const { return type; }
    int  getOrderPara() const { return orderPara; }
    void  setOrderPara(int newOrder) { orderPara=newOrder; }


};

#endif // LIPID_H
