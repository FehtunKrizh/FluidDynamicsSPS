#ifndef CELL_H
#define CELL_H

#include "Material.h"

enum geometry{NoGeometry,Plunger,Die,Powder,Foil,Thermocouple,Powder_1,Powder_2};

//class Border
//{
//    public:
//        bool leftX;
//        bool leftY;
//        bool leftZ;

//        bool rightX;
//        bool rightY;
//        bool rightZ;
//        Border(bool inlx=false, bool inly=false, bool inlz=false,
//               bool inrx=false, bool inry=false, bool inrz=false)
//            :leftX(inlx),leftY(inly),leftZ(inlz),rightX(inrx),rightY(inry),rightZ(inrz){}
//};

template <typename value>
class Cell
{
    public:
        bool Borders;
        double S;
        double rp;
        geometry Geometry;
        value Value;
        Materials Material;
        //Cell(bool inBorder=false):Border(inBorder){}
};


class PhysValue
{
    public:
        double T;
        double Phi;
};



#endif // CELL_H
