/* main.cpp */

//#include <stdio.h>
#include <math.h>
//#include <stdlib.h>
#include <array>
#include <fstream>
#include <iostream>
#include <tuple>
#include <cmath>
//#include <string> //allerede included i iostream
#include <list>
using namespace std;

//Data structre
//struct Position {
//    double x;
//    double y;
//    double z;
//}; //}pos; kan man også skrive? hvad er pos så?
//
//struct fVariables {
//    double x;
//    double y;
//    double z;
//    double vx;
//    double vy;
//    double vz;
//}; //}pos; kan man også skrive? hvad er pos så?

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

     double (*fx)(double,double,double,double,double,double,double,double,struct fParamaters);
     double (*fy)(double,double,double,double,double,double,double,double,struct fParamaters);
     double (*fz)(double,double,double,double,double,double,double,double,struct fParamaters);
     double (*fvx)(double,double,double,double,double,double,double,double,struct fParamaters);
     double (*fvy)(double,double,double,double,double,double,double,double,struct fParamaters);
     double (*fvz)(double,double,double,double,double,double,double,double,struct fParamaters);

     double (*Fv)(double, double, double);
 };


//Baseball functions
//Lige nu er BÅDE variables OG functions smame name...
//måske ikke så godt
double fx(double,double,double,double,double,double,double,double,struct fParamaters);
double fx(double x,double y,double z,double vx,double vy,double vz,double t,double Fv,struct fParamaters fparam){
    double fx;
    fx = vx;
    return fx;
}
double fy(double,double,double,double,double,double,double,double,struct fParamaters);
double fy(double x,double y,double z,double vx,double vy,double vz,double t,double Fv,struct fParamaters fparam){
    double fy;
    fy = vy;
    return fy;
}
double fz(double,double,double,double,double,double,double,double,struct fParamaters);
double fz(double x,double y,double z,double vx,double vy,double vz,double t,double Fv,struct fParamaters fparam){
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

double fvx(double,double,double,double,double,double,double,double,struct fParamaters);
double fvx(double x,double y,double z,double vx,double vy,double vz,double t,double Fv,struct fParamaters fparam){
    double afvx;
    double B = fparam.B;
    double w = fparam.w;
    double phi = fparam.phi;
    double v = sqrt(vx*vx+vy*vy+vz*vz);
    afvx = -Fv*v*vx + B*w*(vz*sin(phi)-vy*cos(phi));
    return afvx;
}
double fvy(double,double,double,double,double,double,double,double,struct fParamaters);
double fvy(double x,double y,double z,double vx,double vy,double vz,double t,double Fv,struct fParamaters fparam){
    double B = fparam.B;
    double w = fparam.w;
    double phi = fparam.phi;
    double v = sqrt(vx*vx+vy*vy+vz*vz);
    double afvy;
    afvy = -Fv*v*vy + B*w*vx*cos(phi);
    return afvy;
}
double fvz(double,double,double,double,double,double,double,double,struct fParamaters);
double fvz(double x,double y,double z,double vx,double vy,double vz,double t,double Fv,struct fParamaters fparam){
    double B = fparam.B;
    double w = fparam.w;
    double phi = fparam.phi;
    double g = fparam.g;
    double v = sqrt(vx*vx+vy*vy+vz*vz);
    double afvz;
    afvz = -g-Fv*v*vz - B*w*vx*sin(phi);
    return afvz;
}




void func1(){cout << "func1 here \n";}
void func2(){cout << "func2 here \n";}
//Baseball RungeKutta
std::tuple<double, double,
            double, double,
            double, double>
RungeKutta2(double,double,double,double,
           double,double,double,double,
           struct fParamaters, functionDaemon);

std::tuple
<double,double,double,double,double,double>
RungeKutta2(
double x,
double y,
double z,
double vx,
double vy,
double vz,
double t,
double h,
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

    double (*fx)(double,double,double,double,double,double,double,double,struct fParamaters);
    double (*fy)(double,double,double,double,double,double,double,double,struct fParamaters);
    double (*fz)(double,double,double,double,double,double,double,double,struct fParamaters);
    double (*fvx)(double,double,double,double,double,double,double,double,struct fParamaters);
    double (*fvy)(double,double,double,double,double,double,double,double,struct fParamaters);
    double (*fvz)(double,double,double,double,double,double,double,double,struct fParamaters);
    fx = f.fx;
    fy = f.fy;
    fz = f.fz;
    fvx = f.fvx;
    fvy = f.fvy;
    fvz = f.fvz;


    //Jeg bør lave en array
    //k1array
    //k2array
    //k3array
    //k4array
    const int neqs = 6;
    double k1a[neqs] = {};
    double k2a[neqs] = {};
    double k3a[neqs] = {};
    double k4a[neqs] = {};
    double k5rk[neqs] = {};

    //Mathematically, alle eqs skal tage samme antal arguments,
    //så det fixer vi lige..
    double (*nfs[neqs])(double,double,double,double,double,double,double,double,struct fParamaters) = {};

    //ftest(*);
    nfs[0] = fx;
    nfs[1] = fy;
    nfs[2] = fz;
    nfs[3] = fvx;
    nfs[4] = fvy;
    nfs[5] = fvz;


    //Doing "krok"
//    k1a[0] = h*fx(x,y,z,vx,vy,vz,t,Fvcalc,fparam);
//    k1a[1] = h*fy(x,y,z,vx,vy,vz,t,Fvcalc,fparam);
//    k1a[2] = h*fz(x,y,z,vx,vy,vz,t,Fvcalc,fparam);
//    k1a[3] = h*fvx(x,y,z,vx,vy,vz,t,Fvcalc,fparam);
//    k1a[4] = h*fvy(x,y,z,vx,vy,vz,t,Fvcalc,fparam);
//    k1a[5] = h*fvz(x,y,z,vx,vy,vz,t,Fvcalc,fparam);

    for (int n=0;n<neqs;++n){
        k1a[n] = h*nfs[n](x,y,z,vx,vy,vz,t,Fvcalc,fparam);
    }

    k1x = h*fx(x,y,z,vx,vy,vz,t,Fvcalc,fparam);
    k1y = h*fy(x,y,z,vx,vy,vz,t,Fvcalc,fparam);
    k1z = h*fz(x,y,z,vx,vy,vz,t,Fvcalc,fparam);
    k1vx = h*fvx(x,y,z,vx,vy,vz,t,Fvcalc,fparam);
    k1vy = h*fvy(x,y,z,vx,vy,vz,t,Fvcalc,fparam);
    k1vz = h*fvz(x,y,z,vx,vy,vz,t,Fvcalc,fparam);

    v = sqrt((vx+0.5*k1vx)*(vx+0.5*k1vx) + (vy+0.5*k1vy)*(vy+0.5*k1vy) + (vz+0.5*k1vz)*(vz+0.5*k1vz));
    Fvcalc = Fv(v,vd,inc);

    k2x = h*fx(x,y,z,vx+0.5*k1vx,vy+0.5*k1vy,vz+0.5*k1vz,t,Fvcalc,fparam);
    k2y = h*fy(x,y,z,vx+0.5*k1vx,vy+0.5*k1vy,vz+0.5*k1vz,t,Fvcalc,fparam);
    k2z = h*fz(x,y,z,vx+0.5*k1vx,vy+0.5*k1vy,vz+0.5*k1vz,t,Fvcalc,fparam);
    k2vx = h*fvx(x,y,z,vx+0.5*k1vx,vy+0.5*k1vy,vz+0.5*k1vz,t,Fvcalc,fparam);
    k2vy = h*fvy(x,y,z,vx+0.5*k1vx,vy+0.5*k1vy,vz+0.5*k1vz,t,Fvcalc,fparam);
    k2vz = h*fvz(x,y,z,vx+0.5*k1vx,vy+0.5*k1vy,vz+0.5*k1vz,t,Fvcalc,fparam);

//    k2a[0] = h*fx(x,y,z,vx+0.5*k1vx,vy+0.5*k1vy,vz+0.5*k1vz,t,Fvcalc,fparam);
//    k2a[1] = h*fy(x,y,z,vx+0.5*k1vx,vy+0.5*k1vy,vz+0.5*k1vz,t,Fvcalc,fparam);
//    k2a[2] = h*fz(x,y,z,vx+0.5*k1vx,vy+0.5*k1vy,vz+0.5*k1vz,t,Fvcalc,fparam);
//    k2a[3] = h*fvx(x,y,z,vx+0.5*k1vx,vy+0.5*k1vy,vz+0.5*k1vz,t,Fvcalc,fparam);
//    k2a[4] = h*fvy(x,y,z,vx+0.5*k1vx,vy+0.5*k1vy,vz+0.5*k1vz,t,Fvcalc,fparam);
//    k2a[5] = h*fvz(x,y,z,vx+0.5*k1vx,vy+0.5*k1vy,vz+0.5*k1vz,t,Fvcalc,fparam);

    for (int n=0;n<neqs;++n){
        k2a[n] = h*nfs[n](x+0.5*k1a[0],y+0.5*k1a[1],z+0.5*k1a[2],
                        vx+0.5*k1a[3],vy+0.5*k1a[4],vz+0.5*k1a[5],t,Fvcalc,fparam);
    }
    //3rd "krok"
    v = sqrt((vx+0.5*k2vx)*(vx+0.5*k2vx) + (vy+0.5*k2vy)*(vy+0.5*k2vy) + (vz+0.5*k2vz)*(vz+0.5*k2vz));
    Fvcalc = Fv(v,vd,inc);

    k3x = h*fx(x,y,z,vx+0.5*k2vx,vy+0.5*k2vy,vz+0.5*k2vz,t,Fvcalc,fparam);
    k3y = h*fy(x,y,z,vx+0.5*k2vx,vy+0.5*k2vy,vz+0.5*k2vz,t,Fvcalc,fparam);
    k3z = h*fz(x,y,z,vx+0.5*k2vx,vy+0.5*k2vy,vz+0.5*k2vz,t,Fvcalc,fparam);
    k3vx = h*fvx(x,y,z,vx+0.5*k2vx,vy+0.5*k2vy,vz+0.5*k2vz,t,Fvcalc,fparam);
    k3vy = h*fvy(x,y,z,vx+0.5*k2vx,vy+0.5*k2vy,vz+0.5*k2vz,t,Fvcalc,fparam);
    k3vz = h*fvz(x,y,z,vx+0.5*k2vx,vy+0.5*k2vy,vz+0.5*k2vz,t,Fvcalc,fparam);

   for (int n=0;n<neqs;++n){
        k3a[n] = h*nfs[n](x+0.5*k2a[0],y+0.5*k2a[1],z+0.5*k2a[2],
                        vx+0.5*k2a[3],vy+0.5*k2a[4],vz+0.5*k2a[5],t,Fvcalc,fparam);
   }

    //4rd "krok"
    v = sqrt((vx+k3vx)*(vx+k3vx) + (vy+k3vy)*(vy+k3vy) + (vz+k3vz)*(vz+k3vz));
    Fvcalc = Fv(v,vd,inc);

    k4x = h*fx(x,y,z,vx+k3vx,vy+k3vy,vz+k3vz,t,Fvcalc,fparam);
    k4y = h*fy(x,y,z,vx+k3vx,vy+k3vy,vz+k3vz,t,Fvcalc,fparam);
    k4z = h*fz(x,y,z,vx+k3vx,vy+k3vy,vz+k3vz,t,Fvcalc,fparam);
    k4vx = h*fvx(x,y,z,vx+k3vx,vy+k3vy,vz+k3vz,t,Fvcalc,fparam);
    k4vy = h*fvy(x,y,z,vx+k3vx,vy+k3vy,vz+k3vz,t,Fvcalc,fparam);
    k4vz = h*fvz(x,y,z,vx+k3vx,vy+k3vy,vz+k3vz,t,Fvcalc,fparam);

   for (int n=0;n<neqs;++n){
        k4a[n] = h*nfs[n](x+k3a[0],y+k3a[1],z+k3a[2],
                        vx+k3a[3],vy+k3a[4],vz+k3a[5],t,Fvcalc,fparam);
   }


    //cout << f(x,y) << " ";
//    x1rk = x + fakt*(k1x+2*k2x+2*k3x+k4x);
//    y1rk = y + fakt*(k1y+2*k2y+2*k3y+k4y);
//    z1rk = z + fakt*(k1z+2*k2z+2*k3z+k4z);
//    vx1rk = vx + fakt*(k1vx+2*k2vx+2*k3vx+k4vx);
//    vy1rk = vy + fakt*(k1vy+2*k2vy+2*k3vy+k4vy);
//    vz1rk = vz + fakt*(k1vz+2*k2vz+2*k3vz+k4vz);

    //cout << f(x,y) << " ";
//    x1rk = x + fakt*(k1a[0]+2*k2a[0]+2*k3a[0]+k4a[0]);
//    y1rk = y + fakt*(k1a[1]+2*k2a[1]+2*k3a[1]+k4a[1]);
//    z1rk = z + fakt*(k1a[2]+2*k2a[2]+2*k3a[2]+k4a[2]);
//    vx1rk = vx + fakt*(k1a[3]+2*k2a[3]+2*k3a[3]+k4a[3]);
//    vy1rk = vy + fakt*(k1a[4]+2*k2a[4]+2*k3a[4]+k4a[4]);
//    vz1rk = vz + fakt*(k1a[5]+2*k2a[5]+2*k3a[5]+k4a[5]);

    double k5rkhelper[neqs] = {x,y,z,vx,vy,vz};

    for (int n=0;n<neqs;++n){
        k5rk[n] = k5rkhelper[n] + fakt*(k1a[n]+2*k2a[n]+2*k3a[n]+k4a[n]);
   }
    //cout << k1 << " " << k2 << "\n" ;
    //x1 = 4.2;
    //y1 = 5.4;
    //cout << x1 << " " << y1 << "\n";

    //return std::make_tuple(x1rk, y1rk, z1rk, vx1rk, vy1rk, vz1rk);
    return std::make_tuple(k5rk[0], k5rk[1], k5rk[2], k5rk[3], k5rk[4], k5rk[5]);
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



struct Foo
{
int id;
};

//void getResult(Foo** fooPtrArray)
int getResult(Foo* fooPtrArray[], unsigned n)
{
cout << "I am in getResult" << endl;
//Foo* fooPtrArray1[2];
//fooPtrArray1 = fooPtrArray;
}

void (*pf[])(void) = {func1, func2};
void test(int jump_index)
{
    /* Call the function specified by jump_index */
    pf[jump_index]();
}

int main(){
    //Foo* fooPtrArray[2] = ;
    //getResult(fooPtrArray, 2);


    test(1);

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
    double k = 1;
    double dt = 0.00001;
    double dt2;
    double eps = 0.047;
    double dx,dy,dz,dvx,dvy,dvz;
    dx=dy=dz=dvx=dvy=dvz=0.01;
    double dxtest,dytest,dztest,dvxtest,dvytest,dvztest;

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

    //tworze parameters do functions
    struct fParamaters fparam;
    fparam.phi = phi;
    fparam.vd = vd;
    fparam.inc = inc;
    fparam.w = w;
    fparam.g = g;
    fparam.B = B;

    //19 w tym chwila?
    int DynamicChanges = 0;

    const int Adaptive = 1;
    switch(Adaptive){
    case 1:

        cout << "Adaptive Runge Kutta \n";
        for (int n=0;n<nmax;++n){
            tie(x2a,y2a,z2a,vx2a,vy2a,vz2a) = RungeKutta2(x[n],y[n],z[n],vx[n],vy[n],vz[n],t[n],dt,fparam,f); //fx,fy,fz,fvx,fvy,fvz
            tie(x1a,y1a,z1a,vx1a,vy1a,vz1a) = RungeKutta2(x[n],y[n],z[n],vx[n],vy[n],vz[n],t[n],0.5*dt,fparam,f);

            dx = std::abs(x2a-x1a);
            dy = std::abs(y2a-y1a);
            dz = std::abs(z2a-z1a);
            dvx = std::abs(vx2a-vx1a);
            dvy = std::abs(vy2a-vy1a);
            dvz = std::abs(vz2a-vz1a);

            //cout << "dy" << dy;
            //cout << "y2 =" << y2a << " "
            //     << "y1 = " << y1a << "\n";

            //cout << "x2 =" << x2a << " "
            //     << "x1 = " << x1a << "\n";


            //albo (dx <<
            //w tym chwila, dx zawiesc..
            //mosliwy robic loop
            //AllTrue=True
            //for dxi in errorlist:
            //if dxi < eps:
            //Pass
            //else:
            //AllTrue=False
            //Break
            if (dx < eps && dy < eps && dz < eps && dvx < eps && dvy < eps && dvz < eps){
            //error is smaller than eps, we accept the step forward
                //cout << "Error was smaller than eps \n";
                //hvis dy < eps, skal vi da også lige prøve
                //at se om vi kan gøre dt STØRRE.
                //Men, måske først EFTER dette step er lavet, så?

                while (dx < eps && dy < eps && dz < eps && dvx < eps && dvy < eps && dvz < eps) {
                    //skal måske lige lave en ny dt2 til testing
                    //så jeg ikke ALTID overwrite den gamle dt...
                    //ellers så overwrite vi den måske eng ang for meget
                    dt2 = 2.0*dt;

                    tie(x2a,y2a,z2a,vx2a,vy2a,vz2a) = RungeKutta2(x[n],y[n],z[n],vx[n],vy[n],vz[n],t[n],dt2,fparam,f);
                    tie(x1a,y1a,z1a,vx1a,vy1a,vz1a) = RungeKutta2(x[n],y[n],z[n],vx[n],vy[n],vz[n],t[n],dt,fparam,f);

                    dxtest = std::abs(x2a-x1a);
                    dytest = std::abs(y2a-y1a);
                    dztest = std::abs(z2a-z1a);
                    dvxtest = std::abs(vx2a-vx1a);
                    dvytest = std::abs(vy2a-vy1a);
                    dvztest = std::abs(vz2a-vz1a);
                    //cout << "doubled dt = " << dt2
                    //     << " dytest = " << dytest << "\n";
                    if (dxtest < eps && dytest < eps && dztest < eps && dvxtest < eps && dvytest < eps && dvztest < eps){
                        //dy jeszcze wieksze niz eps po wyrosla dt
                        dx = dxtest;
                        dz = dztest;
                        dvx = dvxtest;
                        dvy = dvytest;
                        dvz = dvztest;

                        dy = dytest; //nowy dy
                        dt = dt2; //nowy dt
                        //cout << "Made dt bigger";
                        //cout << "\n";
                        //cout << dt;
                        //cout << dytest;
                        //break;
                        DynamicChanges++;
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
                //blad byl mniejszy niz eps, wiec uzywam x1,y1,etc
                //x1a,y1a, etc pochodzi od while loop.
                //nie musisz obliczac x1a,y1a snowu/jeszcze raz
                t[n+1] = t[n]+dt;
                x[n+1] = x1a;
                y[n+1] = y1a;
                z[n+1] = z1a;
                vx[n+1] = vx1a;
                vy[n+1] = vy1a;
                vz[n+1] = vz1a;
            }

            //jesli jeden z nich jest wiekze niz eps, trzeba robic dt mniejszy
            else if (dx > eps || dy > eps || dz > eps || dvx > eps || dvy > eps || dvz > eps)
            //error is bigger than eps, we must make dt smaller
            {
                //dt = 0.5*dt;
                //cout << "Made dt smaller";
                //cout << "\n";

                //cout << "Error was larger than eps \n";
                while (dx > eps || dy > eps || dz > eps || dvx > eps || dvy > eps || dvz > eps){
                    dt2 = 0.5*dt;
                    tie(x2a,y2a,z2a,vx2a,vy2a,vz2a) = RungeKutta2(x[n],y[n],z[n],vx[n],vy[n],vz[n],t[n],dt,fparam,f);
                    tie(x1a,y1a,z1a,vx1a,vy1a,vz1a) = RungeKutta2(x[n],y[n],z[n],vx[n],vy[n],vz[n],t[n],dt2,fparam,f);


                    dxtest = std::abs(x2a-x1a);
                    dztest = std::abs(z2a-z1a);
                    dvxtest = std::abs(vx2a-vx1a);
                    dvytest = std::abs(vy2a-vy1a);
                    dvztest = std::abs(vz2a-vz1a);
                    dytest = std::abs(y2a-y1a);
                    //cout << "Made dt smaller";
                    //cout << "\n";
                    if (dxtest > eps || dytest > eps || dztest > eps || dvxtest > eps || dvytest > eps || dvztest > eps){
                        //The errors(Or just one of them) is still larger than eps, so we accept
                        //the new smaller timestep, and now we have new errors, so we do dx = dxtest
                        dx = dxtest;
                        dz = dztest;
                        dvx = dvxtest;
                        dvy = dvytest;
                        dvz = dvztest;

                        dy = dytest; //this should break the while loop
                        dt = dt2;
                        //cout << "Made dt smaller";
                        //cout << "\n";
                        DynamicChanges++;
                        }
                    else{
                        dt = dt2;
                        break;
                    }
                }

                //we are done decreasing dt, now we take the final step
                //this is perhaps not needed, maybe we can just reuse the old
                //x2a,y2a from before? Or maybe it is needed
                tie(x2a,y2a,z2a,vx2a,vy2a,vz2a) = RungeKutta2(x[n],y[n],z[n],vx[n],vy[n],vz[n],t[n],dt,fparam,f);
                tie(x1a,y1a,z1a,vx1a,vy1a,vz1a) = RungeKutta2(x[n],y[n],z[n],vx[n],vy[n],vz[n],t[n],0.5*dt,fparam,f);

                dx = std::abs(x2a-x1a);
                dy = std::abs(y2a-y1a);
                dz = std::abs(z2a-z1a);
                dvx = std::abs(vx2a-vx1a);
                dvy = std::abs(vy2a-vy1a);
                dvz = std::abs(vz2a-vz1a);
                if (dx < eps && dy < eps && dz < eps && dvx < eps && dvy < eps && dvz < eps)
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
            cout << n << " " << x[n] << " " << y[n] << "\n";
            //cout << x[n+1] << " " << exps[n+1] << " " << y[n+1] << "\n";
            //cout << x[n+1] << " " << v[n+1] << "\n";

        }

//    cout << x[0] << " " << y[0] << "\n";
//    cout << x[1] << " " << y[1] << "\n";
//    cout << x[2] << " " << y[2] << "\n";
//    cout << x[3] << " " << y[3] << "\n";
    cout << "Changes to dt: " << DynamicChanges;
    break;
    case 0:
        cout << "Static Runge Kutta \n";
        dt = dt*10*10*2;
        for (int n=0;n<nmax;++n){
            //struct Position pos;
            //pos.x = x[n];
            //pos.y = y[n];
            //pos.z = z[n];

            tie(x2a,y2a,z2a,vx2a,vy2a,vz2a) = RungeKutta2(x[n],y[n],z[n],vx[n],vy[n],vz[n],t[n],dt,fparam,f);

            t[n+1] = t[n]+dt;
            x[n+1] = x2a;
            y[n+1] = y2a;
            z[n+1] = z2a;
            vx[n+1] = vx2a;
            vy[n+1] = vy2a;
            vz[n+1] = vz2a;
            cout << n << " " << x[n] << " " << y[n] << "\n";
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
