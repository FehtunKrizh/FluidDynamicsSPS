#ifndef MATERIAL_H
#define MATERIAL_H
#include <math.h>


///все параметры представлены в СИ
/// теплоёмкость Дж/(кг*К)
/// плотность кг/м^3
/// теплопроводность Вт/(м*К)
/// теплоёмкость на плотность Дж/(К*м^3)
enum typeMaterial{NoMaterial,Copper,Steel,Alumina,Graphite,Fe,Cu,Ag,Au,Ni,Co};

class Materials
{
    public:
        double k;//теплопроводность
        double cp;//теплоёмкость
        double rho;//плотность
        double sigma;//удельная электропроводность
        typeMaterial type;
        Materials(double ink=1.0,long  double incp=1.0,long  double inrho=1.0,long  double insigma=1.0)
            :k(ink),cp(incp),rho(inrho),sigma(insigma){}
        void refreshMaterial(double T)
        {
            if(this->type==Copper)
            {
                copper(T);
            }
            if(this->type==Alumina)
            {
                alumina(T);
            }
            if(this->type==Steel)
            {
                steel(T);
            }
            if(this->type==Graphite)
            {
                graphite(T);
            }
            if(this->type==Fe)
            {
                mFe(T);
            }
            if(this->type==Cu)
            {
                mCu(T);
            }
            if(this->type==Ag)
            {
                mAg(T);
            }
            if(this->type==Au)
            {
                mAu(T);
            }
            if(this->type==Ni)
            {
                mNi(T);
            }
            if(this->type==Co)
            {
                mCo(T);
            }
            if(this->type==NoMaterial)
            {
                noMaterial(T);
            }
        }

        void alumina(double T)
        {
            //Z.A. Munir Fundamental investigation on the spark plasma sintering synthesis process II Modeling of current and temperature
            double a=7770.25;
            double b=249.4;
            double c=790.15;
            double d=249;
            double e=0.008;

            this->cp=a*T/(b+T)+c*T/(d+T)+e*T;

            a=65181330.4;
            b=-669628.8;
            c=8175.85;
            this->k=(a+T)/(b+c*T);
            this->sigma=1e-8;
            this->rho=3950;//wiki
        }
        void copper(double T)
        {
            //Z.A. Munir Fundamental investigation on the spark plasma sintering synthesis process II Modeling of current and temperature
            this->k=420.66+0.07*T;
            this->cp=355.3+0.1*T;
            this->sigma=1.e9/(5.5+0.038*T);
            //Медь МСр1 ГОСТ 16130-90
            this->rho=8900.0;//кг/м^3
        }
        void graphite(double T)
        {
            this->k=65-0.017*T;
            this->cp=310.5+1.09*T;
            //Eugene A. Olevsky Fundamental Aspects of Spark Plasma Sintering: II. Finite Element Analysis of Scalability
            this->sigma=3.0/(6.083*1e-6+1.4585*1e-6*T/1000+2.3568*1e-6/(T/1000));
            //Eugene A. Olevsky et al Fundamental Aspects of Spark Plasma Sintering II Finit Elemet Analysus of Scalability
            this->rho=1850.0;//кг/м^3
        }
        void steel(double T)
        {
            //Z.A. Munir Fundamental investigation on the spark plasma sintering synthesis process II Modeling of current and temperature
            this->k=11.36+0.0136*T;
            this->cp=446.5+0.162*T;
            this->sigma=1.0/(5.17*1.e-7+6.8*1.e-10*T);
            //С235-С375 ГОСТ 27772-88 google
            this->rho=7850.0;//кг/м^3
        }
        void mFe(double T)
        {
            double d=8.103-0.23761*exp((T-278.23101)/1024.42744);//г/см^3
            this->rho=d*1000;//кг/м^3
            if(T<1042)
            {
                double sqrT=T*T;
                double sqrSQRT=sqrT*sqrT;
                double inter,b1,b2,b3,b4,b5,b6,b7;
                inter=-276.42;
                b1=7.87375;
                b2=-0.03898;
                b3=1.12512e-4;
                b4=-1.95325e-7;
                b5=2.04623e-10;
                b6=-1.199976e-13;
                b7=3.05205e-17;

                if(T<1000.00001)
                {
                    this->cp=inter+
                            b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT+
                            b5*T*sqrSQRT+b6*sqrT*sqrSQRT+
                            b7*T*sqrT*sqrSQRT;
                }
                else
                {
                    this->cp=-9355.85238+10.33095*T;
                }

                inter=877.47694;
                b1=-9.54497;
                b2=0.04857;
                b3=-1.35856e-4;
                b4=2.23425e-7;
                b5=-2.15858e-10;
                b6=1.13449e-13;
                b7=-2.50333e-17;
                this->k=inter+
                        b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT+
                        b5*T*sqrSQRT+b6*sqrT*sqrSQRT+
                        b7*T*sqrT*sqrSQRT;


                inter=-62.30229;
                b1=0.89157;
                b2=-0.00485;
                b3=1.41931e-5;
                b4=-2.28933e-8;
                b5=2.06573e-11;
                b6=-9.61579e-15;
                b7=1.7688e-18;
                double r=inter+
                        b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT+
                        b5*T*sqrSQRT+b6*sqrT*sqrSQRT+
                        b7*T*sqrT*sqrSQRT;
                r=r*1E-8;
                this->sigma=1/r;


            }
            else
            {
                if(T<1183)
                {
                    double sqrT=T*T;
                    double sqrSQRT=sqrT*sqrT;
                    double inter,b1,b2,b3,b4,b5,b6,b7;

                    double y0,A1,t1;
                    y0=710.7273;
                    A1=2.54872e18;
                    t1=29.07;
                    this->cp=A1*exp(-T/t1)+y0;

                    double R0;
                    y0=30.0039;
                    A1=-2.3091e23;
                    R0=-0.05016;
                    this->k=y0+A1*exp(R0*T);

                    inter=-62.30229;
                    b1=0.89157;
                    b2=-0.00485;
                    b3=1.41931e-5;
                    b4=-2.28933e-8;
                    b5=2.06573e-11;
                    b6=-9.61579e-15;
                    b7=1.7688e-18;
                    double r=inter+
                            b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT+
                            b5*T*sqrSQRT+b6*sqrT*sqrSQRT+
                            b7*T*sqrT*sqrSQRT;
                    r=r*1E-8;
                    this->sigma=1/r;

                }
                else
                {
                    if(T<1200)
                    {
                        double y0,A1,t1,x0;
                        y0=607.1;
                        x0=1182.999;
                        A1=109.1;
                        t1=0.65753;
                        this->cp=y0+A1*exp(-(T-x0)/t1);

                        double A2,p,LOGx0;
                        A1=29;
                        A2=34;
                        LOGx0=1400;
                        p=0.07859;
                        this->k=A1+(A2-A1)/(1+pow(10,(LOGx0-T)*p));

                        double A,P,xc;
                        y0=109;
                        xc=1183.001;
                        A=0.7131;
                        P=0.45177;
                        double r=y0+A*pow(abs(T-xc),P);
                        r=r*1E-8;
                        this->sigma=1/r;
                    }
                    else
                    {
                        if(T<1667)
                        {
                            double y0,A1,t1,x0;
                            y0=-1893.548;
                            x0=1200;
                            A1=2500.64;
                            t1=15204.10211;
                            this->cp=y0+A1*exp((T-x0)/t1);

                            double A2,p,LOGx0;
                            A1=29;
                            A2=34;
                            LOGx0=1400;
                            p=0.07859;
                            this->k=A1+(A2-A1)/(1+pow(10,(LOGx0-T)*p));

                            double A,P,xc;
                            y0=109;
                            xc=1183.001;
                            A=0.7131;
                            P=0.45177;
                            double r=y0+A*pow(abs(T-xc),P);
                            r=r*1E-8;
                            this->sigma=1/r;
                        }
                        else
                        {
                            if(T<1811)
                            {
                                this->cp=-230.47664+0.57221*T;

                                this->k=35;

                                double r=61.34+0.0379*T;
                                r=r*1E-8;
                                this->sigma=1/r;
                            }
                            else
                            {
                                this->cp=833.94179+5.29103e-4*T;
                                this->k=39;

                                double r=85.08+0.02646*T;
                                r=r*1E-8;
                                this->sigma=1/r;
                            }
                        }
                    }
                }
            }
        }
        void mCu(double T)
        {
            double sqrT=T*T;
            double sqrSQRT=sqrT*sqrT;
            double inter,b1,b2,b3,b4,b5,b6,b7;
            double y0,A1,t1,x0;
            double r,d;

            if(T<1357.6)
            {

                inter=393.402;
                b1=-0.50309;
                b2=0.0033;
                b3=-8.5995e-6;
                b4=1.2222e-8;
                b5=-9.78889e-12;
                b6=4.1463e-15;
                b7=-7.14965e-19;

                this->cp=inter+
                        b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT+
                        b5*T*sqrSQRT+b6*sqrT*sqrSQRT+
                        b7*T*sqrT*sqrSQRT;
                d=9.09481-5.35635e-4*T;//гр/см^3
                this->rho=d*1000;//кг/м^3


                inter=709.73001;
                b1=-3.68495;
                b2=0.01818;
                b3=-4.77421e-5;
                b4=7.11146e-8;
                b5=-6.04178e-11;
                b6=2.726e-14;
                b7=-5.06466e-18;

                this->k=inter+
                        b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT+
                        b5*T*sqrSQRT+b6*sqrT*sqrSQRT+
                        b7*T*sqrT*sqrSQRT;

                inter=-0.22153;
                b1=0.00475;
                b2=1.46384e-5;
                b3=-4.87661e-8;
                b4=8.48501e-11;
                b5=-7.88318e-14;
                b6=3.77011e-17;
                b7=-7.264e-21;
                r=inter+
                        b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT+
                        b5*T*sqrSQRT+b6*sqrT*sqrSQRT+
                        b7*T*sqrT*sqrSQRT;
                r=r*1e-8;
                this->sigma=1/r;


            }
            else
            {
                this->cp=513.9-5.577e-16*T;

                y0=7.98;
                x0=1357.599;
                A1=0.381;
                t1=3.17823;
                d=y0+A1*exp(-(T-x0)/t1);
                this->rho=d*1000;

                y0=184;
                x0=1357.599;
                A1=133;
                t1=0.46327;
                this->k=y0+A1*exp(-(T-x0)/t1);

                r=7.50832+0.00994*T;
                r=r*1e-8;
                this->sigma=1/r;
            }
        }
        void mAg(double T)
        {
            double d,r;
            double sqrT=T*T;
            double sqrSQRT=sqrT*sqrT;
            double inter,b1,b2,b3,b4,b5,b6,b7,b8;
            if(T<1235)
            {
                this->cp=225.81743+0.02259*T+2.88305e-5*T*T-7.90685e-10*T*T*T;
                inter=9.84694;
                b1=0.01064;
                b2=-6.63267e-5;
                b3=2.32285e-7;
                b4=-5.30449e-10;
                b5=7.77636e-13;
                b6=-6.90137e-16;
                b7=3.33052e-19;
                b8=-6.65903e-23;
                d=inter+
                        b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT+
                        b5*T*sqrSQRT+b6*sqrT*sqrSQRT+
                        b7*T*sqrT*sqrSQRT+b8*sqrSQRT*sqrSQRT;
                this->rho=d*1000;

                inter=417.05409;
                b1=0.12972;
                b2=-3.85634e-4;
                b3=3.18325e-7;
                b4=-1.00119e-10;
                this->k=inter+
                        b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT;
                r=-0.44842+0.00693*T;
                r=r*1e-8;
                this->sigma=1/r;

            }
            else
            {
                this->cp=310.2-2.60335e-15*T;
                d=9.69445-3.031777e-4*T;
                this->rho=d*1000;
                this->k=113.00729+0.03825*T;
                r=6.87453+0.00844*T;
                r=r*1e-8;
                this->sigma=1/r;
            }
        }
        void mAu(double T)
        {
            double d,r;
            double sqrT=T*T;
            double sqrSQRT=sqrT*sqrT;
            double inter,b1,b2,b3,b4,b5,b6,b7,b8;
            inter=217.78864;
            b1=1.51608;
            b2=-0.00832;
            b3=2.277e-5;
            b4=-3.57513e-8;
            b5=3.31764e-11;
            b6=-1.79008e-14;
            b7=5.18246e-18;
            b8=-6.23303e-22;
            this->k=inter+
                    b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT+
                    b5*T*sqrSQRT+b6*sqrT*sqrSQRT+
                    b7*T*sqrT*sqrSQRT+b8*sqrSQRT*sqrSQRT;
            if(T>1400)
            {
                this->k=249.09094;
            }
            if(T<1337.58)
            {
                inter=131.607;
                b1=-0.04812;
                b2=1.82176e-4;
                b3=-2.00277e-7;
                b4=7.96361e-11;
                this->cp=inter+
                        b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT;
                d=19.62876-9.95815e-4*T;
                this->rho=d*1000;
                r=-0.78588+0.01008*T;
                r=r*1e-8;
                this->sigma=1/r;
            }
            else
            {
                this->cp=147.18457+22.31543*exp(-(T-1337.5801)/319.50882);
                d=17.38429-1.60205e-4*T;
                this->rho=d*1000;
                r=11.96172+0.01429*T;
                r=r*1e-8;
                this->sigma=1/r;
            }

        }
        void mNi(double T)
        {
            double d,r;
            double sqrT=T*T;
            double sqrSQRT=sqrT*sqrT;
            double inter,b1,b2,b3,b4,b5,b6;
            if(T<1727.999)
            {
                inter=9.40187;
                b1=-0.00331;
                b2=8.02778e-6;
                b3=-1.01311e-8;
                b4=5.87591e-12;
                b5=-1.2769e-15;
                d=inter+
                        b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT+
                        b5*T*sqrSQRT;
                this->rho=d*1000;
            }
            else
            {
                d=7.85;
                this->rho=d*1000;
            }

            if(T<629.6)
            {
                inter=-18.50325;
                b1=1.3724;
                b2=-0.00568;
                b3=9.4135e-6;
                b4=-5.62787e-9;
                this->k=inter+
                        b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT;
            }
            else
            {
                if(T<1727.999)
                {
                    inter=115.35442;
                    b1=-0.29187;
                    b2=4.85007e-4;
                    b3=-3.0307e-7;
                    b4=6.5396e-11;
                    this->k=inter+
                            b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT;
                }
                else
                {
                    this->k=69;
                }
            }
            if(T<500)
            {
                this->cp=323.1+0.402*T;
            }
            else
            {
                if(T<629.6)
                {
                    this->cp=514.18023+2.2933e-4*exp(0.02133*T);
                }
                else
                {
                    if(T<800)
                    {
                        this->cp=528.91799+1.27221e14*exp(-0.04372*T);
                    }
                    else
                    {
                        if(T<1728)
                        {
                            inter=5567.2987;
                            b1=-24.02819;
                            b2=0.04625;
                            b3=-4.64317e-5;
                            b4=2.5988e-8;
                            b5=-7.74069e-12;
                            b6=9.62246e-16;
                            this->cp=inter+
                                    b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT+
                                    b5*T*sqrSQRT+b6*sqrT*sqrSQRT;
                        }
                        else
                        {
                            this->cp=735-5.6842e-16*T;
                        }
                    }
                }
            }
            if(T<1728.001)
            {
                inter=5.16812;
                b1=-0.06893;
                b2=3.66059e-4;
                b3=-4.50164e-7;
                b4=2.35256e-10;
                b5=-4.51081e-14;
                r=inter+
                        b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT+
                        b5*T*sqrSQRT;
                r=r*1e-8;
                this->sigma=1/r;
            }
            else
            {
                if(T<1728.0011)
                {
                    r=-2.07359e7-12000*T;
                    r=r*1e-8;
                    this->sigma=1/r;
                }
                else
                {
                    if(T<1800)
                    {
                        r=791.01024-0.40973*T;
                        r=r*1e-8;
                        this->sigma=1/r;
                    }
                    else
                    {
                        if(T<2000)
                        {
                            r=-258.8+0.1735*T;
                            r=r*1e-8;
                            this->sigma=1/r;
                        }
                        else
                        {
                            r=88.2;
                            r=r*1e-8;
                            this->sigma=1/r;
                        }
                    }
                }
            }
        }
        void mCo(double T)
        {
            double r,d;
            double sqrT=T*T;
            double sqrSQRT=sqrT*sqrT;
            double inter,b1,b2,b3,b4,b5;
            d=9.06753-6.97402e-4*T;
            this->rho=d*1000;
            if(T<1394)
            {
                inter=-38.01113;
                b1=3.33259;
                b2=-0.00919;
                b3=1.30065e-5;
                b4=-8.86689e-9;
                b5=2.37789e-12;
                this->cp=inter+
                            b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT+
                            b5*T*sqrSQRT;
            }
            else
            {
                if(T<1766.999)
                {
                    double y0=641.93348;
                    double x0=1394;
                    double A1=378.06652;
                    double t1=43.144;
                    this->cp=y0+A1*exp(-(T-x0)/t1);
                }
                else
                {
                    this->cp=688-5.49166e-16*T;
                }
            }
            if(T<1200)
            {
                inter=142.59374;
                b1=-0.20484;
                b2=1.74829e-4;
                b3=-8.91932e-8;
                b4=2.73202e-11;
                this->k=inter+
                            b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT;
            }
            else
            {
                if(T<1767)
                {
                    this->k=300.08185-0.3498*T+1.18097e-4*T*T;
                }
                else
                {
                    this->k=43;
                }
            }
            if(T<1400)
            {
                inter=-1.58566;
                b1=0.02426;
                b2=-1.11624e-5;
                b3=6.3771e-8;
                b4=-2.55989e-11;
                r=inter+
                        b1*T+b2*sqrT+b3*T*sqrT+b4*sqrSQRT;
                r=r*1e-8;
                this->sigma=1/r;
            }
            else
            {
                if(T<1600)
                {
                    r=243.2-0.1115*T;
                    r=r*1e-8;
                    this->sigma=1/r;
                }
                else
                {
                    if(T<1767)
                    {
                        r=64.8+3.9537e-91*exp(0.11999*T);
                        r=r*1e-8;
                        this->sigma=1/r;
                    }
                    else
                    {
                        r=67.03586+0.02549*T;
                        r=r*1e-8;
                        this->sigma=1/r;
                    }
                }
            }
        }

        void noMaterial(double T)
        {
            this->k=-0.00001+0.0*T;
            this->cp=0.00001+0.0*T;
            this->sigma=-0.00001+0.0*T;
            this->rho=0.000001;
        }
};





#endif // MATERIAL_H
