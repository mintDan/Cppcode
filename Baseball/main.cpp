/* iteration1.c */

//#include <stdio.h>
#include <math.h>
//#include <stdlib.h>
#include <array>
#include <fstream>
#include <iostream>
#include <tuple>
#include <cmath>
#include <string>
#include <list>
using namespace std;


double stuff(double);
double stuff(double a)
{
    printf("Hi from stuff");
    printf("a = %3d",a);
}

void ptest(int, double *);
void ptest(int n, double *fact)
{
    *fact = 1;
    printf("n = %3d",n);

}





double f1(double, double,double);
double f1(double t, double x, double y)
{
    double fx;
    fx = y;
    //cout << "f here" << fxy;
    //cout << fxy;
    return fx;

}


double f2(double,double,double);
double f2(double t, double x, double y)
{
    double fv,k;
    k = 1;
    fv = -k*x;
    //cout << "f here" << fxy;
    //cout << fxy;
    return fv;

}

//lidt usikker på arguments tilo runge kutta, men will find out
double RungeKutta(double,double,double,double);
double RungeKutta(double t, double x, double y, double h)
{
    double k1,k2,k3,k4,y1;
    double fakt = 0.1666666;
    k1 = h*f1(t,x,y);
    k2 = h*f1(t,x+h/2,y+k1/2);
    k3 = h*f1(t,x+h/2,y+k2/2);
    k4 = h*f1(t,x+h,y+k3);
    //cout << f(x,y) << " ";
    y1 = y + fakt*(k1+2*k2+2*k3+k4);
    //cout << k1 << " " << k2 << "\n" ;
    //cout << y1 << "\n";
    return y1;
}


std::tuple<double, double> RungeKutta1(double,double,double,double);
std::tuple<double, double> RungeKutta1(double t, double x, double y, double h)
{
    //jeg har her x+h/2....
    //MEN... h er jo dt... og x er position, så units passer ikke!
    //Det er noget med der er en k1x,k1y, right?
    //Hvis jeg husker rigtigt...
    //Så hmmm... måske lidt som partial derivatives
    //først step man x+k1x, og holder y constant?
    //så step man y+k1y og holder x constant? etc maybe?
    double k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,x1,y1;
    double fakt = 0.1666666;

    k1x = h*f1(t,x,y);
    k1y = h*f2(t,x,y);

    k2x = h*f1(t+h/2,x+0.5*k1x,y+0.5*k1y);
    k2y = h*f2(t+h/2,x+0.5*k1x,y+0.5*k1y);

    k3x = h*f1(t+h/2,x+0.5*k2x,y+0.5*k2y);
    k3y = h*f2(t+h/2,x+0.5*k2x,y+0.5*k2y);

    k4x = h*f1(t+h,x+k3x,y+k3y);
    k4y = h*f2(t+h,x+k3x,y+k3y);
    //cout << f(x,y) << " ";
    x1 = x + fakt*(k1x+2*k2x+2*k3x+k4x);
    y1 = y + fakt*(k1y+2*k2y+2*k3y+k4y);
    //cout << k1 << " " << k2 << "\n" ;
    //x1 = 4.2;
    //y1 = 5.4;
    //cout << x1 << " " << y1 << "\n";

    return std::make_tuple(x1, y1);
}


//Baseball functions
//Lige nu er BÅDE variables OG functions smame name...
//måske ikke så godt
double fx(double);
double fx(double vx){
    double fx;
    fx = vx;
    return fx;
}
double fy(double);
double fy(double vy){
    double fy;
    fy = vy;
    return fy;
}
double fz(double);
double fz(double vz){
    double fz;
    fz = vz;
    return fz;
}


double Fv(double, double, double);
double Fv(double v,double vd, double inc){
    double Fv;
    Fv = 0.0039 + 0.0058/(1+exp((v-vd)/inc));
    return Fv;
}

double fvx(double,double,double,double,double,double,double,double);
double fvx(double Fv,double v,double vx,double B, double w, double vz,double phi,double vy){
    double afvx;
    afvx = -Fv*v*vx + B*w*(vz*sin(phi)-vy*cos(phi));
    return afvx;
}
double fvy(double, double, double, double, double,double,double,double);
double fvy(double Fv,double v, double vx,double B, double w, double vz,double phi,double vy){
    double afvy;
    afvy = -Fv*v*vy + B*w*vx*cos(phi);
    return afvy;
}
double fvz(double,double,double,double,double,double,double,double,double);
double fvz(double Fv,double v, double vx,double B, double w, double vz,double phi,double vy,double g){
    double afvz;
    afvz = -g-Fv*v*vz - B*w*vx*sin(phi);
    return afvz;
}



//Data structre
struct Position {
    double x;
    double y;
    double z;
}pos;


struct fParamaters {
    double phi;
    double vd;
    double inc;
    double w;
    double g;
    double B;
}lol;

// struct of function pointers
//https://stackoverflow.com/questions/18663521/function-pointer-in-struct
struct functionDaemon{
     int id;
     //double (*funcp); // function pointer
     //double  (*fp)(double);      // Function pointer

     double (*fx)(double);
     double (*fy)(double);
     double (*fz)(double);
     double (*fvx)(double,double,double,double,double,double,double,double);
     double (*fvy)(double,double,double,double,double,double,double,double);
     double (*fvz)(double,double,double,double,double,double,double,double,double);

     double (*Fv)(double, double, double);
 };


//Baseball RungeKutta
std::tuple<double, double,
            double, double,
            double, double>
RungeKutta2(double,double,double,double,
           double,double,double,double,struct fParamaters, functionDaemon);

std::tuple
<double,double,
double,double,
double,double>
RungeKutta2(double t,
             double x,
              double y,
            double h,
            double z,
             double vx,
             double vy,
              double vz,
//              double (*f1)(double),
//            double (*f2)(double),
//            double (*f3)(double),
//            double (*f4)(double,double,double,double,double,double,double,double),
//            double (*f5)(double,double,double,double,double,double,double,double),
//            double (*f6)(double,double,double,double,double,double,double,double,double),
            fParamaters fparam,
            functionDaemon f){

    double k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,x1rk,y1rk,z1rk;
    double vx1rk,vy1rk,vz1rk;
    double k1z,k2z,k3z,k4z;
    double k1vx,k2vx,k3vx,k4vx;
    double k1vy,k2vy,k3vy,k4vy;
    double k1vz,k2vz,k3vz,k4vz;
    double Fvcalc;
    double v;
    double fakt = 0.1666666;

    double phi = fparam.phi;
    double vd = fparam.vd;
    double inc = fparam.inc;
    double w = fparam.w;
    double g = fparam.g;
    double B = fparam.B;

    v = sqrt(vx*vx + vy*vy + vz*vz);
    Fvcalc = Fv(v,vd,inc);

    //Jeg bør have de fleste functions
    //fra struct, også den sidste Fv fra pointer,
    //således jeg har code example

    double (*f1)(double);
    double (*f2)(double);
    double (*f3)(double);
    double (*fvx)(double,double,double,double,double,double,double,double);
    double (*fvy)(double,double,double,double,double,double,double,double);
    double (*fvz)(double,double,double,double,double,double,double,double,double);
    f1 = f.fx;
    f2 = f.fy;
    f3 = f.fz;
    fvx = f.fvx;
    fvy = f.fvy;
    fvz = f.fvz;

    //Doing steps
    k1x = h*f1(vx);
    k1y = h*f2(vy);
    k1z = h*f3(vz);
    k1vx = h*fvx(Fvcalc,v,vx,B,w,vz,phi,vy);
    k1vy = h*fvy(Fvcalc,v,vx,B,w,vz,phi,vy);
    k1vz = h*fvz(Fvcalc,v,vx,B,w,vz,phi,vy,g);

    v = sqrt((vx+0.5*k1vx)*(vx+0.5*k1vx) + (vy+0.5*k1vy)*(vy+0.5*k1vy) + (vz+0.5*k1vz)*(vz+0.5*k1vz));
    Fvcalc = Fv(v,vd,inc);

    k2x = h*fx(vx+0.5*k1vx);
    k2y = h*fy(vy+0.5*k1vy);
    k2z = h*fz(vz+0.5*k1vz);
    k2vx = h*fvx(Fvcalc,v,vx+0.5*k1vx,B,w,vz+0.5*k1vz,phi,vy+0.5*k1vy);
    k2vy = h*fvy(Fvcalc,v,vx+0.5*k1vx,B,w,vz+0.5*k1vz,phi,vy+0.5*k1vy);
    k2vz = h*fvz(Fvcalc,v,vx+0.5*k1vx,B,w,vz+0.5*k1vz,phi,vy+0.5*k1vy,g);

    v = sqrt((vx+0.5*k2vx)*(vx+0.5*k2vx) + (vy+0.5*k2vy)*(vy+0.5*k2vy) + (vz+0.5*k2vz)*(vz+0.5*k2vz));
    Fvcalc = Fv(v,vd,inc);

    k3x = h*fx(vx+0.5*k2vx);
    k3y = h*fy(vy+0.5*k2vy);
    k3z = h*fz(vz+0.5*k2vz);
    k3vx = h*fvx(Fvcalc,v,vx+0.5*k2vx,B,w,vz+0.5*k2vz,phi,vy+0.5*k2vy);
    k3vy = h*fvy(Fvcalc,v,vx+0.5*k2vx,B,w,vz+0.5*k2vz,phi,vy+0.5*k2vy);
    k3vz = h*fvz(Fvcalc,v,vx+0.5*k2vx,B,w,vz+0.5*k2vz,phi,vy+0.5*k2vy,g);

    v = sqrt((vx+k3vx)*(vx+k3vx) + (vy+k3vy)*(vy+k3vy) + (vz+k3vz)*(vz+k3vz));
    Fvcalc = Fv(v,vd,inc);

    k4x = h*fx(vx+k3vx);
    k4y = h*fy(vy+k3vy);
    k4z = h*fz(vz+k3vz);
    k4vx = h*fvx(Fvcalc,v,vx+k1vx,B,w,vz+k1vz,phi,vy+k1vy);
    k4vy = h*fvy(Fvcalc,v,vx+k1vx,B,w,vz+k1vz,phi,vy+k1vy);
    k4vz = h*fvz(Fvcalc,v,vx+k1vx,B,w,vz+k1vz,phi,vy+k1vy,g);

    //cout << f(x,y) << " ";
    x1rk = x + fakt*(k1x+2*k2x+2*k3x+k4x);
    y1rk = y + fakt*(k1y+2*k2y+2*k3y+k4y);
    z1rk = z + fakt*(k1z+2*k2z+2*k3z+k4z);
    vx1rk = vx + fakt*(k1vx+2*k2vx+2*k3vx+k4vx);
    vy1rk = vy + fakt*(k1vy+2*k2vy+2*k3vy+k4vy);
    vz1rk = vz + fakt*(k1vz+2*k2vz+2*k3vz+k4vz);
    //cout << k1 << " " << k2 << "\n" ;
    //x1 = 4.2;
    //y1 = 5.4;
    //cout << x1 << " " << y1 << "\n";

    return std::make_tuple(x1rk, y1rk, z1rk, vx1rk, vy1rk, vz1rk);
}
//const int nmax = 200;


/* Set max. allowable no. of iterations */
//#define NITER 30

//det skal faktisk ikke være double, men arrays
void PrintToFile(
double *t,
double *x,
double *y,
double *z,
double *vx,
double *vy,
double *vz,
int nsize,
string filestring){
        //det skal være HELE arrays som jeg pass til PrintToFile
    ofstream myfile;
    myfile.open(filestring);
    for (int n=1;n<nsize;++n){
        myfile << t[n] << " ";
        myfile << x[n] << " ";
        myfile << y[n] << " ";
        myfile << z[n] << " ";
        myfile << vx[n] << " ";
        myfile << vy[n] << " ";
        myfile << vz[n] << " " << "\n";
    }
    myfile.close();
}


void func1(){cout << "func1 here \n";}
void func2(){cout << "func2 here \n";}

//typedef void (*lolpointer[2])(); works kinda halfway
//void(*lolpointer[2])(); works kinda halfway

//http://forums.codeguru.com/showthread.php?456631-list-of-function-pointers
//try this method now
typedef void (*fp)();

//Denne her function, skal vel bare tage en list?
//no?
//void ftest(void *lolpointer[2]){
//    //lolpointer[0]();
//    //lolpointer[1]();
//    //(*lolpointer[0])();
//}

void ftest(void *lolpointer){
    //cout << (*lolpointer[0])();
    lolpointer;
    cout << "\n";
    cout << lolpointer;
    &lolpointer;
    cout << "\n";
    cout << &lolpointer;
    //&lolpointer[0];

    //(*lolpointer[0])();
    //(*lolpointer)[0]();
    //(*lolpointer[0])();
}

void ftest2(int *lolpointer){
    cout << lolpointer[0] << "\n";
    cout << lolpointer[1] << "\n";
    cout << lolpointer[2] << "\n";
}



int main(){

    struct functionDaemon f;
    f.fx = fx;
    f.fy = fy;
    f.fz = fz;
    f.fvx = fvx;
    f.fvy = fvy;
    f.fvz = fvz;
    f.Fv = Fv;
    //const lolpointer ar[2] = {&func1, &func2};
    //ponizej dziala tez z const lolpointer
    //lolpointer[0] = func1;
     //lolpointer[1] = func2;
     //lolpointer[0]();
     //lolpointer[1]();
    //ar[0]();
    //ar[1]();
    //ftest((*lolpointer));
    int heythere[3] = {88,99,122};
    //int heythere2[2] = {func1,*func2}; dette virker ikke
    //Så, man kan ikke bare lave en simple list af function pointers
    ftest2(heythere);
    //BEGGE DE HER TO VIRKER!!!
    //Så endnu en forskel her.
    //void (*lolpointer[2])() = {&func1,&func2};
    void (*lolpointer[2])() = {func1,func2};

    //ftest(*);
    lolpointer[0]();
    lolpointer[1]();
    ftest(lolpointer);
    //ahh, <fp> er nok "type" af list?
    //og fpl er name af list, prolly?
    //problemet er, at vi sender det jo ikke ind i en function..
    std::list<fp> fpl;
	fpl.push_back(&func1); //http://www.cplusplus.com/reference/vector/vector/push_back/ push_back append til list?
	fpl.push_back(&func2);

    cout << "\n";
    cout << reinterpret_cast<void*>(ftest);
    cout << "\n";
    double k = 1,dt = 0.00001;
    double dt2;
    double dytest;
    double eps = 0.001;
    double dy = 0.01;
    double v0 = 37.9984;
    double theta = 0.0174532925;
    double phi = 0;
    double vd = 35;
    double inc = 5;
    double w = 188.495559;
    double g = 9.81;
    double B = 4.1/(10*10*10*10);
    double L = 18.44;
    double T = L/v0;

//    array<double,nmax> x {0};
//    array<double,nmax> v {sqrt(k)};
//    array<double,nmax> t {0};
//    array<double,nmax> y {0};
//    array<double,nmax> z {0};
//    array<double,nmax> vx {v0*cos(theta)};
//    array<double,nmax> vy {0};
//    array<double,nmax> vz {v0*sin(theta)};
    const int nmax = 250;//T/dt;
    double t[nmax] = {0};
    double x[nmax] = {0};
    double v[nmax] = {sqrt(k)};
    double y[nmax] = {0};
    double z[nmax] = {0};
    double vx[nmax] = {v0*cos(theta)};
    double vy[nmax] = {0};
    double vz[nmax] = {v0*sin(theta)};

    double x1a,y1a,z1a,x2a,y2a,z2a;
    double vx1a,vy1a,vz1a,vx2a,vy2a,vz2a;

    struct fParamaters fparam;
    fparam.phi = phi;
    fparam.vd = vd;
    fparam.inc = inc;
    fparam.w = w;
    fparam.g = g;
    fparam.B = B;

    const int Adaptive = 0;
    switch(Adaptive){
    case 1:
        cout << "Adaptive Runge Kutta \n";
        for (int n=0;n<nmax;++n){
            //struct Position pos;
            //pos.x = x[n];
            //pos.y = y[n];
            //pos.z = z[n];
            //cout << "timestep n = " << n << "\n";

            tie(x2a,y2a,z2a,vx2a,vy2a,vz2a) = RungeKutta2(t[n],x[n],y[n],dt,z[n],vx[n],vy[n],vz[n],fparam,f); //fx,fy,fz,fvx,fvy,fvz
            tie(x1a,y1a,z1a,vx1a,vy1a,vz1a) = RungeKutta2(t[n],x[n],y[n],0.5*dt,z[n],vx[n],vy[n],vz[n],fparam,f);


            dy = std::abs(y2a-y1a);
            //cout << "dy" << dy;
            //cout << "y2 =" << y2a << " "
            //     << "y1 = " << y1a << "\n";

            //cout << "x2 =" << x2a << " "
            //     << "x1 = " << x1a << "\n";
            if (dy < eps){
            //error is smaller than eps, we accept the step forward
                //cout << "Error was smaller than eps \n";
                //hvis dy < eps, skal vi da også lige prøve
                //at se om vi kan gøre dt STØRRE.
                //Men, måske først EFTER dette step er lavet, så?

                while (dy < eps) {
                    //skal måske lige lave en ny dt2 til testing
                    //så jeg ikke ALTID overwrite den gamle dt...
                    //ellers så overwrite vi den måske eng ang for meget
                    dt2 = 2.0*dt;

                    tie(x2a,y2a,z2a,vx2a,vy2a,vz2a) = RungeKutta2(t[n], x[n], y[n],dt2,z[n],vx[n], vy[n], vz[n],fparam,f);
                    tie(x1a,y1a,z1a,vx1a,vy1a,vz1a) = RungeKutta2(t[n], x[n], y[n],dt,z[n],vx[n], vy[n], vz[n],fparam,f);

                    dytest = std::abs(y2a-y1a);

                    //cout << "doubled dt = " << dt2
                    //     << " dytest = " << dytest << "\n";
                    if (dytest < eps){
                        //dy jeszcze wieksze niz eps po wyrosla dt
                        dy = dytest; //nowy dy
                        dt = dt2; //nowy dt
                        //cout << "Made dt bigger";
                        //cout << "\n";
                        //cout << dt;
                        //cout << dytest;
                        //break;
                    }
                    else{
                        //ostatni dytest byl wieksze niz eps, wiec
                        //zatrzymam loop z break
                        //ponadto, nie przydzielam dt = dt2.
                        //Dlatego, powinno dzialac dt.
                        //dt = dt2*0.5;

                        //cout << "error no longer smaller than eps \n";
                        break;
                    }
                }
                //blad byl mniejszy niz eps, wiec uzywam x2,y2,etc
                t[n+1] = t[n]+dt;
                x[n+1] = x1a;
                y[n+1] = y1a;
                z[n+1] = z1a;
                vx[n+1] = vx1a;
                vy[n+1] = vy1a;
                vz[n+1] = vz1a;
            }


            else if (dy > eps)
            //error is bigger than eps, we must make dt smaller
            {
                //dt = 0.5*dt;
                //cout << "Made dt smaller";
                //cout << "\n";

                //cout << "Error was larger than eps \n";
                while (dy>eps){
                    dt2 = 0.5*dt;
                    tie(x2a,y2a,z2a,vx2a,vy2a,vz2a) = RungeKutta2(t[n], x[n], y[n],dt,z[n],vx[n], vy[n], vz[n],fparam,f);
                    tie(x1a,y1a,z1a,vx1a,vy1a,vz1a) = RungeKutta2(t[n], x[n], y[n],dt2,z[n],vx[n], vy[n], vz[n],fparam,f);

                    dytest = std::abs(y2a-y1a);
                    //cout << "Made dt smaller";
                    //cout << "\n";
                    if (dytest > eps){
                        dy = dytest; //this should break the while loop
                        dt = dt2;
                        //cout << "Made dt smaller";
                        //cout << "\n";
                        }
                    else{
                        dt = dt2;
                        break;
                    }
                }
                tie(x2a,y2a,z2a,vx2a,vy2a,vz2a) = RungeKutta2(t[n], x[n], y[n],dt,z[n],vx[n], vy[n], vz[n],fparam,f);
                tie(x1a,y1a,z1a,vx1a,vy1a,vz1a) = RungeKutta2(t[n], x[n], y[n],dt*0.5,z[n],vx[n], vy[n], vz[n],fparam,f);

                dy = std::abs(y2a-y1a);
                if (dy < eps)
                //now we should have corrected for dy, it should now be lower than eps
                {
                    t[n+1] = t[n]+dt;
                    x[n+1] = x2a;
                    y[n+1] = y2a;
                    z[n+1] = z2a;
                    vx[n+1] = vx2a;
                    vy[n+1] = vy2a;
                    vz[n+1] = vz2a;
                }
            }

            //exps[n+1] = exp(x[n]);
            cout << x[n] << " " << y[n] << "\n";
            //cout << x[n+1] << " " << exps[n+1] << " " << y[n+1] << "\n";
            //cout << x[n+1] << " " << v[n+1] << "\n";
        }

//    cout << x[0] << " " << y[0] << "\n";
//    cout << x[1] << " " << y[1] << "\n";
//    cout << x[2] << " " << y[2] << "\n";
//    cout << x[3] << " " << y[3] << "\n";
    break;
    case 0:
        cout << "Static Runge Kutta \n";
        dt = dt*10*10*2;
        for (int n=0;n<nmax;++n){
            //struct Position pos;
            //pos.x = x[n];
            //pos.y = y[n];
            //pos.z = z[n];

            tie(x2a,y2a,z2a,vx2a,vy2a,vz2a) = RungeKutta2(t[n],x[n],y[n],dt,z[n],vx[n],vy[n],vz[n],fparam,f);

            t[n+1] = t[n]+dt;
            x[n+1] = x2a;
            y[n+1] = y2a;
            z[n+1] = z2a;
            vx[n+1] = vx2a;
            vy[n+1] = vy2a;
            vz[n+1] = vz2a;
            cout << x[n] << " " << y[n] << "\n";
        }
    break;
    }

    //napisam do plik
    int nsize = nmax;
    string filestring = "Slider.txt";
    //cout << filestring;
    PrintToFile(t,x,y,z,vx,vy,vz,nsize,filestring);

    return 0;
}
