#include <time.h>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <malloc.h>
#include <iostream>

#include "Cell.h"
#include "Material.h"

using namespace std;

const double gamma=0.2;//параметр сходимости для элиптического уравнения(та часть которая решает потенциалл)
const unsigned int SizeX=32+1; //количество ячеек по x;
const unsigned int SizeY=32; //количество ячеек по y;
const unsigned int SizeZ=32+2; //количество ячеек по z;
const unsigned int Size=SizeX*SizeY*SizeZ;
const unsigned int SizeIndex=(SizeX-2)*(SizeY)*(SizeZ-2);

unsigned int center[SizeIndex];

unsigned int leftX[SizeIndex];
unsigned int leftY[SizeIndex];
unsigned int leftZ[SizeIndex];
unsigned int rightX[SizeIndex];
unsigned int rightY[SizeIndex];
unsigned int rightZ[SizeIndex];


//линейные размеры области
const double lx=0.02;
const double ly=2*M_PI;
const double lz=0.08;

const double Rd=0.02;//радиус пресс формы
const double Rp=0.01;//радиус плунжеров
const double Hd=0.04;//высота пресс формы
const double Hp=0.02;//длина на которыю выступает наружу плунжеры
const double Hm=0.012;//высота от начала пресс формы до материалла(плунжеры в глубину)

const double dx=lx/(SizeX-1);
const double dy=ly/SizeY;
const double dz=lz/(SizeZ-2);

const double I=500;//I=1000;
const double topU=1.0;
const double bottomU=0.0;

const double dt=0.001;//шаг по времени

double CFLX=dt/(dx*dx);
double CFLY=dt/(dy*dy);
double CFLZ=dt/(dz*dz);


///все параметры представлены в СИ
/// теплоёмкость Дж/(кг*К)
/// плотность кг/м^3
/// теплопроводность Вт/(м*К)
/// теплоёмкость на плотность Дж/(К*м^3)


Cell<PhysValue> Cells1[Size];//на нечетном временном слое
Cell<PhysValue> Cells2[Size];//на четном временном слое


//функия задания индексов(кроме границ)
void initialIndedx(void);
//функция задание геометрии образца
void setupGeometry(Cell<PhysValue> *cells);
//функция задания начальных условий
void initialConditions(Cell<PhysValue> *cells);
//функция сохранения результатов
void saveResults(Cell<PhysValue> *cells, unsigned int time);

//решаем уравнение теплопроводности
void solverHeat(Cell<PhysValue> *cells1, Cell<PhysValue> *cells2);
//решем урвнение теплопроводности без источника тепла
void solverHeatCooling(Cell<PhysValue> *cells1, Cell<PhysValue> *cells2);
//граничные условия для уравнения теплопроводности
void borderHeat(Cell<PhysValue> *cells1, Cell<PhysValue> *cells2);

//решаем уравнение токопроводимости
void solverElectric(Cell<PhysValue> *cells1, Cell<PhysValue> *cells2);
//граничные условия дляуравнения токопроводимости
void borderElectric(Cell<PhysValue> *cells1, Cell<PhysValue> *cells2);

//вспомогательные функции
//расчет потока
double flux(double x,double y);
//потери тепла на границе из-за испускательной способности
double heatLossRadiative(const double T, const double v=1.0);//Вт м^-2
//потери тепла на границе из-за конвекции
double heatLossConvection(const double T, const double k=370.0);
//приток тепла из-за нагрева тока
double heatCurrent(const double current, const double fieldStrength);
//вычисление индекса
unsigned int findIndex(const unsigned int x,const unsigned int y, const unsigned int z);

void loadFile(void);
bool continueCalculation(bool);


int main(void)
{
    FILE *f;
    unsigned int timeOFF=600/dt;
    unsigned int timeColling=400/dt;

//    FILE *h;
//    Cell<PhysValue> a;
//    a.Material.type=Co;
//    h=fopen("Thermo.dat ","w");
//    for(double T=300;T<2000;T=T+0.01)
//    {
//        a.Value.T=T;
//        a.Material.refreshMaterial(T);
//        fprintf(h,"\n%lf\t%lf\t%lf\t%lf\t%lf\t",T,a.Material.cp,a.Material.rho,1e8*1/a.Material.sigma,a.Material.k);
//    }
//    fclose(h);
//    printf("The end");
//    return 0;

    //printf("\nsetupGeometry");
    setupGeometry(Cells1);
    //printf("\ninitialConditions");
    initialConditions(Cells1);

    //printf("\ninitialIndex");
    initialIndedx();


    for(unsigned int i=0;i<Size;++i)
    {
        Cells2[i]=Cells1[i];
    }
    f=fopen("Thermocouple.dat","w");
    fprintf(f,"t\tTp\tTs\tTe\n");
    fclose(f);

    for(unsigned int time=0;time<timeOFF;++time)
    {

        if(time%1000==0)//(time%10000==0)
        {
            printf("\ntime=%lf",time*dt);
            saveResults(Cells2,time);
        }
        if(time<timeColling)
        {
            if(time%2)
            {

                for(unsigned int i=0;i<0.2*sqrt(Size);++i)
                {
                    solverElectric(Cells1,Cells2);
                    borderElectric(Cells1,Cells2);
                }

                solverHeat(Cells1,Cells2);
                borderHeat(Cells1,Cells2);
            }
            else
            {

                for(unsigned int i=0;i<0.2*sqrt(Size);++i)
                {
                    solverElectric(Cells2,Cells1);
                    borderElectric(Cells1,Cells2);
                }

                solverHeat(Cells2,Cells1);
                borderHeat(Cells2,Cells1);

            }
        }
        else
        {
            if(time%2)
            {

                solverHeatCooling(Cells1,Cells2);
                borderHeat(Cells1,Cells2);
            }
            else
            {
                solverHeatCooling(Cells2,Cells1);
                borderHeat(Cells2,Cells1);

            }
        }

    }
    printf("The End.\n");
    return 0;
}

void initialIndedx(void)
{
    unsigned int count=0;
    for(unsigned int x=1;x<SizeX-1;++x)
    {
        for(unsigned int y=0;y<SizeY;++y)
        {
            for(unsigned int z=1;z<SizeZ-1;++z)
            {
                //определение индексов

                center[count]=findIndex(x,y,z);

                leftX[count]=findIndex(x-1,y,z);
                leftY[count]=findIndex(x,y-1,z);
                leftZ[count]=findIndex(x,y,z-1);

                rightX[count]=findIndex(x+1,y,z);
                rightY[count]=findIndex(x,y+1,z);
                rightZ[count]=findIndex(x,y,z+1);

                if(y==0)
                {
                    leftY[count]=findIndex(x,SizeY-1,z);
                }
                if(y==SizeY-1)
                {
                    rightY[count]=findIndex(x,0,z);
                }


                count++;
            }
        }
    }
}

void setupGeometry(Cell<PhysValue> *cells)
{
    for(unsigned int x=0,index;x<SizeX;++x)
    {
        for(unsigned int y=0;y<SizeY;++y)
        {
            for(unsigned int z=0;z<SizeZ;++z)
            {
                index=findIndex(x,y,z);
                cells[index].rp=x*dx;
                cells[findIndex(x,y,z)].Geometry=NoGeometry;

                if((x<SizeX-1)&&(z>0)&&(z<SizeZ-1))
                {
                    if(x*dx<Rp)
                    {
                        cells[index].Geometry=Plunger;
                        cells[index].S=M_PI*Rp*Rp;
                        if((z*dz>Hp)&&((z-1)*dz<=Hp+Hd))
                        {
                            cells[index].Geometry=Plunger;
                            cells[index].S=M_PI*Rd*Rd;
                        }
//                        if(((z-1)*dz>Hp+Hm)&&(z*dz<Hp+Hd-Hm))
//                        {
//                            cells[index].Geometry=Powder;
//                            cells[index].S=M_PI*Rd*Rd;
//                        }
                        if(((z-1)*dz>Hp+Hm)&&(z*dz<Hp+Hd-Hm))
                        {
                            cells[index].Geometry=Powder_1;
                            cells[index].S=M_PI*Rd*Rd;
                        }
                        if(((z-1)*dz>Hp+Hm)&&(z*dz<Hp+Hd-1.5*Hm))
                        {
                            cells[index].Geometry=Powder_2;
                            cells[index].S=M_PI*Rd*Rd;
                        }
                    }
                    if((x*dx>=Rp)&&(x*dx<=Rd))
                    {
                        if((z*dz>Hp)&&((z-1)*dz<=Hp+Hd))
                        {
                            cells[index].Geometry=Die;
                            cells[index].S=M_PI*Rd*Rd;
                        }
                    }
                    if(((x+1)*dx==Rp))
                    {
                        if((z*dz>Hp)&&((z-1)*dz<=Hp+Hd))
                        {
                            cells[index].Geometry=Foil;
                            cells[index].S=M_PI*Rd*Rd;
                        }
                    }
                }
                if(x==SizeX*2/3)
                {
                    if(y==SizeY/2)
                    {
                        if(z==SizeZ/2)
                        {
                            cells[index].Geometry=Thermocouple;
                        }
                    }
                }


//                if((x<SizeX-1)&&(z>0)&&(z<SizeZ-1))
//                {
//                    if(x*dx<0.3*lx)
//                    {
//                        if((z*dz<0.45*lz)||(z*dz>0.55*lz))
//                        {
//                            cells[findIndex(x,y,z)].Geometry=Plunger;
//                        }
//                        else
//                        {
//                            cells[findIndex(x,y,z)].Geometry=Powder;
//                        }
//                    }
//                    else
//                    {
//                        cells[findIndex(x,y,z)].Geometry=Die;
//                    }
//                }

                //определение границ
                cells[index].Borders=false;
                if(x==SizeX-2)
                {
                    cells[index].Borders=true;
                }
                if((z==SizeZ-2)||(z==1))
                {
                    cells[index].Borders=true;
                }
                if((z-1)*dz<=Hp)
                {
                    if((x+1)*dx>=Rp)
                    {
                        cells[index].Borders=true;
                    }
                }
                if(z*dz>=Hp+Hd)
                {
                    if((x+1)*dx>=Rp)
                    {
                        cells[index].Borders=true;
                    }
                }
                if((x+1)*dx>=Rp)
                {
                    if(z*dz<Hp)
                    {
                        cells[index].Borders=true;
                    }
                    if(z*dz>Hp+Hd)
                    {
                        cells[index].Borders=true;
                    }
                }
            }
        }
    }
    //printf("\n asda");


    for(unsigned int i=0;i<Size;++i)
    {
        cells[i].Material.type=NoMaterial;
        if(!(cells[i].Geometry^Die))
        {
            cells[i].Material.type=Graphite;
        }
        if(!(cells[i].Geometry^Foil))
        {
            cells[i].Material.type=Graphite;
        }
        if(!(cells[i].Geometry^Thermocouple))
        {
            cells[i].Material.type=Graphite;
        }
//        if(!(cells[i].Geometry^Powder))
//        {
//            cells[i].Material.type=Alumina;//Steel;
//        }
        if(!(cells[i].Geometry^Powder_1))
        {
            cells[i].Material.type=Alumina;//Copper;
        }
        if(!(cells[i].Geometry^Powder_2))
        {
            cells[i].Material.type=Alumina;
        }
        if(!(cells[i].Geometry^Plunger))
        {
            cells[i].Material.type=Graphite;//Steel;
        }
    }
}


void initialConditions(Cell<PhysValue> *cells)
{
    for(unsigned int x=0,index;x<SizeX;++x)
    {
        for(unsigned int y=0;y<SizeY;++y)
        {
            for(unsigned int z=0;z<SizeZ;++z)
            {
                index=findIndex(x,y,z);
                cells[index].Value.T=cells[index].Value.Phi=0.0;
                cells[index].Value.T=0.0;
                cells[index].Material.refreshMaterial(0.0);

                if(cells[index].Material.type!=NoMaterial)
                {
                    cells[index].Value.T=300.0;
                    cells[index].Material.refreshMaterial(300.0);
                    cells[index].Value.Phi=(z*topU/(SizeZ-1));
                }


                if(cells[index].Geometry==Plunger)
                {
                    if(z==1)
                    {
                        cells[index].Value.Phi=bottomU;//
                    }
                    if(z==SizeZ-2)
                    {
                        cells[index].Value.Phi=topU;
                    }
                }

            }
        }
    }
}

void saveResults(Cell<PhysValue> *cells, unsigned int time)
{
    FILE *outFile;
    char fileName[512];
    sprintf(fileName,"SPS_time=%lf.csv",time*dt);
    outFile=fopen(fileName,"w+");
    if(!outFile)
    {
        printf("File no open outFile, completion of the program");
        return ;
    }
    fprintf(outFile,"x,y,z,T,Phi,Geomatry,Material,Bordes\n");
    for(unsigned int index,x=0;x<SizeX-1;++x)
    {
        for(unsigned int  y=0;y<SizeY;++y)
        {
            for(unsigned int z=1;z<SizeZ-1;++z)
            {
                index=findIndex(x,y,z);
                if(cells[index].Geometry!=NoGeometry)
                {
                    fprintf(outFile,"%lf,%lf,%lf,%lf,%lf,%d,%d,%d\n",
                                                       x*dx,
                                                       y*dy,
                                                       (z-1)*dz,
                                                       cells[index].Value.T,
                                                       cells[index].Value.Phi,
                                                       cells[index].Geometry,
                                                       cells[index].Material.type,
                                                       cells[index].Borders);
                }
            }
        }
    }
    fclose(outFile);
    outFile=fopen("Thermocouple.dat","a");
    double temp=0.0,surfTemp=0.0,therm=0.0;
    double count1=0,count2=0;
    for(unsigned int i=0;i<Size;++i)
    {
        if((cells[i].Geometry==Powder)||(cells[i].Geometry==Powder_1)||(cells[i].Geometry==Powder_2))
        {
            count1++;
            temp+=cells[i].Value.T;
        }
        if(cells[i].Geometry==Die)
        {
            if(cells[i].Borders)
            {
                count2++;
                surfTemp+=cells[i].Value.T;
            }
        }
        if(cells[i].Geometry==Thermocouple)
        {
            therm=cells[i].Value.T;
        }
    }
    fprintf(outFile,"%lf\t%lf\t%lf\t%lf\n",time*dt,temp/count1,surfTemp/count2,therm);
    fclose(outFile);
}


void solverHeat(Cell<PhysValue> *cells1, Cell<PhysValue> *cells2)
{
    double E,J,Ex,Ey,Ez;
    double Q,rRhoCp,rrp,rp;
    double cflX,cflY,cflZ;
    double ke,kw,kn,ks,kt,kb;
    double Tp,Te,Tw,Tn,Ts,Tt,Tb;

    for(unsigned int i=0;i<SizeIndex;++i)
    {

        if((cells2[center[i]].Material.type!=NoMaterial))
        {
            //double ds=(rp[center[i]]*dx*dy);
            J=I/cells2[center[i]].S;

            //<a|b>=g_ija_ib_j
            //длина вектора а |a|=sqrt<a|a>=sqrt(g_ija_ia_j)
            //E=-grad(U)={du/dx,du/dy,du/dz} в цилиндрических {du/dr,1/r du/dphi, du/dz}
            //<E|E>=g_ijE_iE_j g при i=/=j g=0 иначе g(1,r^2,1) <E|E>=(du/dr)^2+r^2*(1/r du/dphi)^2+(du/dz)^2



            Ex=(cells2[rightX[i]].Value.Phi-cells2[leftX[i]].Value.Phi)/(2*dx);
            Ey=(cells2[rightY[i]].Value.Phi-cells2[leftY[i]].Value.Phi)/(2*dy);
            Ez=(cells2[rightZ[i]].Value.Phi-cells2[leftZ[i]].Value.Phi)/(2*dz);

            if(cells2[rightX[i]].Material.type==NoMaterial)
            {
                Ex=(cells2[center[i]].Value.Phi-cells2[leftX[i]].Value.Phi)/dx;
            }
            if(cells2[leftX[i]].Material.type==NoMaterial)
            {
                Ex=(cells2[rightX[i]].Value.Phi-cells2[center[i]].Value.Phi)/dx;
            }

            if(cells2[rightZ[i]].Material.type==NoMaterial)
            {
                Ez=(cells2[center[i]].Value.Phi-cells2[leftZ[i]].Value.Phi)/dz;
            }
            if(cells2[leftZ[i]].Material.type==NoMaterial)
            {
                Ez=(cells2[rightZ[i]].Value.Phi-cells2[center[i]].Value.Phi)/dz;
            }


            E=sqrt(Ex*Ex+Ey*Ey+Ez*Ez);



            Q=heatCurrent(J,E);



            Tp=cells2[center[i]].Value.T;

            cells2[center[i]].Material.refreshMaterial(Tp);


            Tw=cells2[leftX[i]].Value.T;
            Te=cells2[rightX[i]].Value.T;

            kw=flux(cells2[center[i]].Material.k,cells2[leftX[i]].Material.k);
            ke=flux(cells2[center[i]].Material.k,cells2[rightX[i]].Material.k);
            if(cells2[rightX[i]].Material.type==NoMaterial)
            {
                Te=0.0;
                ke=0.0;
            }
            if(cells2[leftX[i]].Material.type==NoMaterial)
            {
                Tw=0.0;
                kw=0.0;
            }



            //тут уже учтенно что индексы centeri=0 следовательно левый будет i=SizeY-1
            Ts=cells2[leftY[i]].Value.T;
            Tn=cells2[rightY[i]].Value.T;

            ks=flux(cells2[center[i]].Material.k,cells2[leftY[i]].Material.k);
            kn=flux(cells2[center[i]].Material.k,cells2[rightY[i]].Material.k);


            Tb=cells2[leftZ[i]].Value.T;
            Tt=cells2[rightZ[i]].Value.T;


            kb=flux(cells2[center[i]].Material.k,cells2[leftZ[i]].Material.k);
            kt=flux(cells2[center[i]].Material.k,cells2[rightZ[i]].Material.k);
            if(cells2[rightZ[i]].Material.type==NoMaterial)
            {
                Tt=0.0;
                kt=0.0;
            }
            if(cells2[leftZ[i]].Material.type==NoMaterial)
            {
                Tb=0.0;
                kb=0.0;
            }

            rp=cells2[center[i]].rp;
            rrp=1.0/(rp*rp);

            rRhoCp=1.0/(cells2[center[i]].Material.rho*cells2[center[i]].Material.cp);
            cflX=CFLX*rRhoCp;
            cflY=CFLY*rRhoCp*rrp;
            cflZ=CFLZ*rRhoCp;



            cells1[center[i]].Value.T=Tp+rRhoCp*Q*dt+cflX*(ke*(Te-Tp)-kw*(Tp-Tw))+cflY*(kn*(Tn-Tp)-ks*(Tp-Ts))+cflZ*(kt*(Tt-Tp)-kb*(Tp-Tb));
            if(cells1[center[i]].Borders)
            {
                if(cells1[center[i]].Geometry==Die)
                {
                    //излучательная способность 0,8 типо матрица окутана шубкой
                    cells1[center[i]].Value.T-=heatLossRadiative(cells1[center[i]].Value.T,0.8)*rRhoCp*dt/(sqrt(dx*dz));
                }
                else
                {
                    if(cells1[center[i]].Geometry==Plunger)
                    {
                        cells1[center[i]].Value.T-=heatLossConvection(cells1[center[i]].Value.T)*rRhoCp*dt/(sqrt(dx*dz));
                    }
                    else
                    {
                        cells1[center[i]].Value.T-=heatLossRadiative(cells1[center[i]].Value.T)*rRhoCp*dt/(sqrt(dx*dz));
                    }
                }
            }

            if(cells1[center[i]].Value.T!=cells1[center[i]].Value.T)
            {
                printf("\n value=%.10lf ",Tp+rRhoCp*Q*dt+cflX*(ke*(Te-Tp)-kw*(Tp-Tw))+cflY*(kn*(Tn-Tp)-ks*(Tp-Ts))+cflZ*(kt*(Tt-Tp)-kb*(Tp-Tb)));
                printf("\nT=%lf",cells1[center[i]].Value.T);
                getchar();
            }
        }
    }
}

void solverHeatCooling(Cell<PhysValue> *cells1, Cell<PhysValue> *cells2)
{
    double rRhoCp,rrp,rp;
    double cflX,cflY,cflZ;
    double ke,kw,kn,ks,kt,kb;
    double Tp,Te,Tw,Tn,Ts,Tt,Tb;

    for(unsigned int i=0;i<SizeIndex;++i)
    {

        if((cells2[center[i]].Material.type!=NoMaterial))
        {

            Tp=cells2[center[i]].Value.T;

            cells2[center[i]].Material.refreshMaterial(Tp);


            Tw=cells2[leftX[i]].Value.T;
            Te=cells2[rightX[i]].Value.T;

            kw=flux(cells2[center[i]].Material.k,cells2[leftX[i]].Material.k);
            ke=flux(cells2[center[i]].Material.k,cells2[rightX[i]].Material.k);
            if(cells2[rightX[i]].Material.type==NoMaterial)
            {
                Te=0.0;
                ke=0.0;
            }
            if(cells2[leftX[i]].Material.type==NoMaterial)
            {
                Tw=0.0;
                kw=0.0;
            }



            //тут уже учтенно что индексы centeri=0 следовательно левый будет i=SizeY-1
            Ts=cells2[leftY[i]].Value.T;
            Tn=cells2[rightY[i]].Value.T;

            ks=flux(cells2[center[i]].Material.k,cells2[leftY[i]].Material.k);
            kn=flux(cells2[center[i]].Material.k,cells2[rightY[i]].Material.k);


            Tb=cells2[leftZ[i]].Value.T;
            Tt=cells2[rightZ[i]].Value.T;


            kb=flux(cells2[center[i]].Material.k,cells2[leftZ[i]].Material.k);
            kt=flux(cells2[center[i]].Material.k,cells2[rightZ[i]].Material.k);
            if(cells2[rightZ[i]].Material.type==NoMaterial)
            {
                Tt=0.0;
                kt=0.0;
            }
            if(cells2[leftZ[i]].Material.type==NoMaterial)
            {
                Tb=0.0;
                kb=0.0;
            }

            rp=cells2[center[i]].rp;
            rrp=1.0/(rp*rp);

            rRhoCp=1.0/(cells2[center[i]].Material.rho*cells2[center[i]].Material.cp);
            cflX=CFLX*rRhoCp;
            cflY=CFLY*rRhoCp*rrp;
            cflZ=CFLZ*rRhoCp;



            cells1[center[i]].Value.T=Tp+cflX*(ke*(Te-Tp)-kw*(Tp-Tw))+cflY*(kn*(Tn-Tp)-ks*(Tp-Ts))+cflZ*(kt*(Tt-Tp)-kb*(Tp-Tb));
            if(cells1[center[i]].Borders)
            {
                if(cells1[center[i]].Geometry==Die)
                {
                    //излучательная способность 0,8 типо матрица окутана шубкой
                    cells1[center[i]].Value.T-=heatLossRadiative(cells1[center[i]].Value.T,0.8)*rRhoCp*dt/(sqrt(dx*dz));
                }
                else
                {
                    if(cells1[center[i]].Geometry==Plunger)
                    {
                        cells1[center[i]].Value.T-=heatLossConvection(cells1[center[i]].Value.T)*rRhoCp*dt/(sqrt(dx*dz));
                    }
                    else
                    {
                        cells1[center[i]].Value.T-=heatLossRadiative(cells1[center[i]].Value.T)*rRhoCp*dt/(sqrt(dx*dz));
                    }
                }
            }
        }
    }
}

void borderHeat(Cell<PhysValue> *cells1, Cell<PhysValue> *cells2)
{
    //r=0;
    for(unsigned int z=1;z<SizeZ-1;++z)
    {
        unsigned int index;
        double T=0.0;
        for(unsigned int y=0;y<SizeY;++y)
        {
            index=findIndex(1,y,z);
            T+=cells2[index].Value.T;
        }
        T/=SizeY;
        for(unsigned int y=0;y<SizeY;++y)
        {
            index=findIndex(0,y,z);
            cells1[index].Value.T=T;
            cells1[index].Material.refreshMaterial(T);
        }
    }
}

void solverElectric(Cell<PhysValue> *cells1, Cell<PhysValue> *cells2)
{
    double rrp,rp;
    double ke,kw,kn,ks,kt,kb;//sigma
    double A,Ue,Uw,Un,Us,Ut,Ub;

    for(unsigned int i=0;i<SizeIndex;++i)
    {
        if(cells2[center[i]].Material.type!=NoMaterial)
        {

            Uw=cells2[leftX[i]].Value.Phi;
            Ue=cells2[rightX[i]].Value.Phi;

            kw=flux(cells2[center[i]].Material.sigma,cells2[leftX[i]].Material.sigma);
            ke=flux(cells2[center[i]].Material.sigma,cells2[rightX[i]].Material.sigma);
            if(cells2[rightX[i]].Material.type==NoMaterial)
            {
                Ue=0.0;
                ke=0.0;
            }
            if(cells2[leftX[i]].Material.type==NoMaterial)
            {
                Uw=0.0;
                kw=0.0;
            }


            //тут уже учтенно что индексы centeri=0 следовательно левый будет i=SizeY-1
            Us=cells2[leftY[i]].Value.Phi;
            Un=cells2[rightY[i]].Value.Phi;

            ks=flux(cells2[center[i]].Material.sigma,cells2[leftY[i]].Material.sigma);
            kn=flux(cells2[center[i]].Material.sigma,cells2[rightY[i]].Material.sigma);


            Ub=cells2[leftZ[i]].Value.Phi;
            Ut=cells2[rightZ[i]].Value.Phi;


            kb=flux(cells2[center[i]].Material.sigma,cells2[leftZ[i]].Material.sigma);
            kt=flux(cells2[center[i]].Material.sigma,cells2[rightZ[i]].Material.sigma);
            if(cells2[rightZ[i]].Material.type==NoMaterial)
            {
                Ut=0.0;
                kt=0.0;
            }
            if(cells2[leftZ[i]].Material.type==NoMaterial)
            {
                Ub=0.0;
                kb=0.0;
            }

            rp=cells2[center[i]].rp;
            rrp=1/(rp*rp);

            A=1.0/((ke+kw)/(dx*dx)+(kn+ks)*rrp/(dy*dy)+(kb+kt)/(dz*dz));


            cells1[center[i]].Value.Phi=A*((ke*Ue+kw*Uw)/(dx*dx)+(kn*Un+ks*Us)*rrp/(dy*dy)+(kb*Ub+kt*Ut)/(dz*dz));
        }
    }
}

void borderElectric(Cell<PhysValue> *cells1, Cell<PhysValue> *cells2)
{

    //r=0;
    for(unsigned int z=1;z<SizeZ-1;++z)
    {
        unsigned int index;
        double U=0.0;
        for(unsigned int y=0;y<SizeY;++y)
        {
            index=findIndex(1,y,z);
            U+=cells2[index].Value.Phi;
        }
        U/=SizeY;
        for(unsigned int y=0;y<SizeY;++y)
        {
            index=findIndex(0,y,z);
            cells1[index].Value.Phi=U;
        }
    }

    for(unsigned int x=0;x<SizeX;++x)
    {
        for(unsigned int y=0,index;y<SizeY;++y)
        {
            index=findIndex(x,y,1);
            if(cells1[index].Geometry==Plunger)
            {
                cells1[index].Value.Phi=bottomU;
            }
            index=findIndex(x,y,SizeZ-2);
            if(cells1[index].Geometry==Plunger)
            {
                cells1[index].Value.Phi=topU;
            }
        }
    }
}


double flux(double x,double y)
{
    if(x<0)
    {
        return 0.0;
    }
    if(y<0)
    {
        return 0.0;
    }
    return 2*x*y/(x+y);//return 0.5*(x+y);
}


double heatLossRadiative(const double T, const double v)
{
    const double sigma=5.670367*1e-8;// постоянная Стефана-Больцмана Вт м^-2 K^-4
    //const double v=1.0;//излучательная способность
    const double T0_4=8100000000;//температура к которой стремится остыть 300
    return sigma*v*(T*T*T*T-T0_4);
}

double heatLossConvection(const double T, const double k)
{
    double T_0=300;
    return k*(T-T_0);
}

double heatCurrent(const double current, const double fieldStrength)
{
    return current*fieldStrength;
}


double sqr(const double x)
{
    return x*x;
}

unsigned int findIndex(const unsigned int x,const unsigned int y, const unsigned int z)
{
    unsigned int rv=x+y*SizeX+z*SizeX*SizeY;
    if(rv>Size)
    {
        return 0;
    }
    else
    {
        return rv;
    }
}
