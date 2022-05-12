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
    int  getOrder() const { return orderPara; }
    void  setOrder(int newOrder) { orderPara=newOrder; }
    void  setID(int newID) { ID=newID; }
    void  setType(int newType) { type=newType; }


};

#endif // LIPID_H
