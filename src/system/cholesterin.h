#ifndef CHOLESTERIN_H
#define CHOLESTERIN_H

class Cholesterin
{
public:
    Cholesterin(int, int ,int,bool);
    int posX, posY, ID;
    bool occupied;
    int  getID()   const { return ID; }
    void  setID(int newID) { ID=newID; }
    
};

#endif // CHOLESTERIN_H
